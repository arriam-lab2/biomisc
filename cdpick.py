#! /usr/bin/env python


"""

A wrapper script around cd-hit-est-2d for closed reference OTU-picking. Performs
global alignment.

    cdpick.py -r ref.fna -o mapping.txt -i samples.fna

"""

import sys

if sys.version_info < (3, 6):
    print('This tool requires Python >= 3.6')
    sys.exit(1)

from typing import TypeVar, Callable, Iterable, Optional, List
from itertools import groupby
import subprocess as sp
import operator as op
import tempfile
import os
import re

import click
from fn import F, _ as var

A = TypeVar('A')
DEVNULL = open(os.devnull, 'w')
ENV = '/usr/bin/env'
QUITE = dict(stdout=DEVNULL, stderr=sp.STDOUT)
CDHIT = 'cd-hit-est-2d'
SEQID = re.compile('>(.+?)\.\.\.').findall


def transform_cluster(noempty, cluster: Iterable[str]) -> Optional[List[str]]:
    seqids = [SEQID(line)[0] for line in cluster]
    return (seqids if len(seqids) > 1 else None) if noempty else seqids


def onpath(executable: str) -> bool:
    """
    Can /usr/bin/env find the executable?
    :param executable: an executable to find
    :return:
    """
    return sp.run([ENV, executable], **QUITE).returncode != 127


def validate(validator: Callable[[A], bool], message: str, ctx, param: str,
             value: A):
    if not validator(value):
        raise click.BadParameter(message, ctx=ctx, param=param)
    return value


@click.command('cdpick', help=__doc__,
               context_settings=dict(help_option_names=['-h', '--help']))
@click.option('-i', '--input', required=True,
              type=click.Path(exists=True, dir_okay=False, resolve_path=True),
              help='A FASTA file to map')
@click.option('-r', '--reference', required=True,
              type=click.Path(exists=True, dir_okay=False, resolve_path=True),
              help='Reference dataset (FASTA or FASTQ)')
@click.option('-o', '--output', required=True,
              type=click.Path(exists=False, dir_okay=False, resolve_path=True),
              callback=F(validate, lambda v: not os.path.exists(v),
                         'output exists'),
              help='Output destination.')
@click.option('-a', '--accurate', is_flag=True, default=False,
              help='Run in accurate mode. By default, a sequence is clustered '
                   'to the first cluster that meet the threshold (fast '
                   'cluster). In this mode the program will cluster it into '
                   'the most similar cluster that meet the threshold. This '
                   'might take several times more runtime, though.')
@click.option('-s', '--similarity', type=float, default=0.97,
              callback=F(validate, lambda v: 0.5 <= v <= 1, 'not in [0, 1]'),
              help='Sequence similarity cutoff value; a floating point number '
                   'within [0.5, 1].')
@click.option('-t', '--threads', type=int, default=1,
              callback=F(validate, var > 0, 'must be positive'),
              help='The number of CPU threads to use.')
@click.option('-m', '--memory', type=int, default=1000,
              callback=F(validate, var >= 100, 'should be at least 100MB'),
              help='Maximum amount of RAM available to CD-HIT (must be at '
                   'least 100MB).')
@click.option('-e', '--supress_empty', is_flag=True, default=False)
def cdpick(input: str, reference: str, output: str, accurate: bool,
           similarity: float, threads: int, memory: int, supress_empty: bool):
    # make sure cdhit is available
    if not onpath(CDHIT):
        raise click.UsageError(
            'No cd-hit-est-2d executable found; is it on your PATH?'
        )
    with tempfile.TemporaryDirectory() as root:
        base = os.path.join(root, os.path.splitext(os.path.basename(output))[0])
        clusters = f'{base}.clstr'
        command = [
            ENV, CDHIT, '-i', reference, '-c', str(similarity),
            '-g', str(int(accurate)), '-T', str(threads), '-M', str(memory),
            '-o', base, '-i2', input
        ]
        process = sp.run(command)
        if process.returncode:
            raise click.ClickException(
                'CD-HIT failed; please, read its error logs for more details'
            )
        with open(clusters) as lines, open(output, 'w') as out:
            clusters = (
                F(map, str.strip) >> (filter, bool) >>
                (lambda x: groupby(x, lambda l: l.startswith('>'))) >>
                (filter, lambda x: not x[0]) >> (map, op.itemgetter(1)) >>
                (map, F(transform_cluster, supress_empty)) >> (filter, bool)
            )(lines)
            for clust in clusters:
                print('\t'.join(clust), file=out)


if __name__ == '__main__':
    cdpick()
