#! /usr/bin/env python


"""

A wrapper script around cd-hit-est-2d for closed reference OTU-picking. Performs
global alignment.

    cdpick.py -r ref.fna -o mapping.txt -i samples.fna

"""

from contextlib import ExitStack
import operator as op
import os
import re
import tempfile
import subprocess as sp
from itertools import groupby, chain
from typing import Iterable, Optional, List, Union, Tuple

from fn import F

from pipeline import util

CDHIT = 'cd-hit-est-2d'
GZIP = 'gzip'
SEQID = re.compile('>(.+?)\.\.\.').findall


# TODO add an import-time warnings about cd-hit's and/or gzip's absence

def transform_cluster(noempty, cluster: Iterable[str]) -> Optional[List[str]]:
    seqids = [SEQID(line)[0] for line in cluster]
    return (seqids if len(seqids) > 1 else None) if noempty else seqids


def pick(reference: str, accurate: bool, similarity: float, threads: int,
         memory: int, supress_empty: bool,
         input: Union[Tuple[str], Tuple[str, str]], output: str):
    # make sure cdhit and gzip are available
    if not util.onpath(CDHIT):
        raise RuntimeError(
            'No cd-hit-est-2d executable found; is it on your PATH?'
        )
    if not util.onpath(GZIP):
        raise RuntimeError(
            'No gzip executable found; is it on your PATH?'
        )
    # check input compression status
    uncompressed = []
    with ExitStack() as stack:
        for path in input:
            if util.isgzipped(path):
                buffer = tempfile.NamedTemporaryFile(dir='')
                stack.enter_context(buffer)
                sp.run([util.ENV, 'gzip', '-df', path], stdout=buffer)
                buffer.flush()
                uncompressed.append(buffer.name)
            else:
                uncompressed.append(path)

        base = os.path.splitext(os.path.basename(output))[0]
        clusters = f'{base}.clstr'
        command = [
            util.ENV, CDHIT, '-i', reference, '-c', str(similarity),
            '-g', str(int(accurate)), '-T', str(threads), '-M', str(memory),
            '-o', base, *chain(*zip(('-i2', '-j2'), uncompressed))
        ]
        process = sp.run(command)
        if process.returncode:
            raise RuntimeError(
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
    raise RuntimeError
