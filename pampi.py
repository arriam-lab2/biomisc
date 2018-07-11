import operator as op
import tempfile
import os
from functools import reduce
from typing import List, Callable, TypeVar

import click
from fn import F, _ as X


TMPDIR = 'tmpdir'
FASTQ = 'FASTQ'
FASTA = 'FASTA'

A = TypeVar('A')


def validate(validator: Callable[[A], bool], message: str, ctx, param: str,
             value: A):
    if not validator(value):
        raise click.BadParameter(message, ctx=ctx, param=param)
    return value


@click.group(chain=True, invoke_without_command=True)
@click.option('-i', '--input', required=True)
@click.option('-t', '--tempdir', default=tempfile.gettempdir(),
              type=click.Path(exists=False, dir_okay=False, resolve_path=True),
              callback=F(validate, lambda v: os.path.exists(v),
                         'tmpdir does not exist'),
              help='Temporary directory location')
def pampi(input, tempdir, pattern, group):
    pass


@pampi.resultcallback()
def pipeline(processors, input, tempdir: str):
    # TODO handle compilation type error
    with tempfile.TemporaryDirectory(dir=tempdir) as root:
        pass


# TODO add validators

@pampi.command('qc')
@click.pass_context
@click.option('-l', '--minlen', type=int)
@click.option('-q', '--minqual', type=str, help='e.g. 18:3')  # TODO document
@click.option('-c', '--head_crop', type=int)
@click.option('-o', '--output', required=True,
              type=click.Path(exists=False, resolve_path=True),
              callback=F(validate, lambda v: not os.path.exists(v),
                         'output destination exists'),
              help='Output destination.')
def qc():
    pass


@pampi.command('join')
@click.pass_context
@click.option('-o', '--output', required=True,
              type=click.Path(exists=False, resolve_path=True),
              callback=F(validate, lambda v: not os.path.exists(v),
                         'output exists'),
              help='Output destination.')
@click.option('-p', '--pattern', type=str,
              help='A Python regular expression. By default the expression '
                   'is used to split basenames and select the first value in '
                   'the split. You might instead want to use a group regex to '
                   'extract the first occurrence of the group by specifying '
                   'the --group flag')
@click.option('--group', is_flag=True, default=False)
def join():
    pass


@pampi.command('pick')
@click.option('-r', '--reference', required=True,
              type=click.Path(exists=True, dir_okay=False, resolve_path=True),
              help='Reference dataset (FASTA or FASTQ)')
@click.option('-a', '--accurate', is_flag=True, default=False,
              help='Run in accurate mode. By default, a sequence is clustered '
                   'to the first cluster that meet the threshold (fast '
                   'cluster). In this mode the program will cluster it into '
                   'the most similar cluster that meet the threshold. This '
                   'might take several times more runtime, though.')
@click.option('-s', '--similarity', type=float, default=0.97,
              callback=F(validate, lambda v: 0.5 <= v <= 1, 'not in [0.5, 1]'),
              help='Sequence similarity cutoff value; a floating point number '
                   'within [0.5, 1].')
@click.option('-t', '--threads', type=int, default=1,
              callback=F(validate, X > 0, 'must be positive'),
              help='The number of CPU threads to use.')
@click.option('-m', '--memory', type=int, default=1000,
              callback=F(validate, X >= 100, 'should be at least 100MB'),
              help='Maximum amount of RAM available to CD-HIT (must be at '
                   'least 100MB).')
@click.option('-e', '--supress_empty', is_flag=True, default=False)
@click.option('-o', '--output', required=True,
              type=click.Path(exists=False, dir_okay=False, resolve_path=True),
              callback=F(validate, lambda v: not os.path.exists(v),
                         'output exists'),
              help='Output destination.')
def pick():
    pass


if __name__ == '__main__':
    pampi()