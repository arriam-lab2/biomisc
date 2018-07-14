import operator as op
import operator as op
import os
import re
import shutil
import subprocess as sp
from itertools import groupby, chain
from typing import Iterable, Optional, List, Union, Tuple

from fn import F

from pipeline.pampi import util, data

CDHIT = 'cd-hit-est-2d'
SEQID = re.compile('>(.+?)\.\.\.').findall
OUT_EXT = '.clstr'


# TODO add an import-time warnings about cd-hit's and/or gzip's absence


def transform_cluster(drop_empty: bool, cluster: Iterable[str]) \
        -> Optional[List[str]]:
    seqids = [SEQID(line)[0] for line in cluster]
    return (seqids if len(seqids) > 1 else None) if drop_empty else seqids


def parse_cdhit_clusters(drop_empty: bool, handle: Iterable[str]) \
        -> Iterable[List[str]]:
        return (
            F(map, str.strip) >> (filter, bool) >>
            (lambda x: groupby(x, lambda l: l.startswith('>'))) >>
            (filter, lambda x: not x[0]) >> (map, op.itemgetter(1)) >>
            (map, F(transform_cluster, drop_empty)) >> (filter, bool)
        )(handle)


# TODO :param supress_empty: suppress empty clusters, i.e. do not report
# TODO                       unobserved references.

@util.fallible(RuntimeError, sp.CalledProcessError)
def cdhit(reference: str, accurate: bool, similarity: float, threads: int,
          memory: int, input: Union[Tuple[str], Tuple[str, str]], output: str) \
        -> Optional[Tuple[str, str]]:
    """
    A basic wrapper around cd-hit-est-2d
    :param reference: path to reference sequences
    :param accurate:  run in accurate mode
    :param similarity: similarity threshold
    :param threads: the number of threads to use
    :param memory: RAM limit; this is a hard limit: cd-hit will stop working
    past that limit.
    :param input: input sequences, either a single FAST[A/Q] file or a pair
    of these (for paired-end libraries); compressed fiiles are not supported.
    :param output: output path for sequences observed in the reference; cd-hit
    will also save a clustering file named '{}.clstr'.format(output).
    :return: path to a fasta file with sequences observed in the reference
    """
    # make sure cdhit is available
    executable = shutil.which(CDHIT)
    if not executable:
        raise RuntimeError(
            'No cd-hit-est-2d executable found; is it on your PATH?'
        )
    clusters = f'{output}.clstr'
    command = [
        executable, '-i', reference, '-c', str(similarity),
        '-g', str(int(accurate)), '-T', str(threads), '-M', str(memory),
        '-o', output, *chain(*zip(['-i2', '-j2'], input))
    ]
    process = sp.run(command)
    if process.returncode:
        raise RuntimeError(
            'CD-HIT failed; please, read its error logs for more details'
        )
    return output, clusters


# TODO add a nondesctructive debug mode?
@util.fallible(RuntimeError, FileNotFoundError)
def cdpick(tmpdir: str, input: data.SampleReads, outdir: Optional[str],
           drop_empty: bool, **cdhit_options) -> Optional[data.SampleClusters]:
    if not os.path.exists(tmpdir):
        raise ValueError(f'temporary directory {tmpdir} does not exist')
    output = (util.randname(tmpdir, OUT_EXT) if outdir is None else
              os.path.join(outdir, input.name+OUT_EXT))
    cdhit_tempout = util.randname(tmpdir, '')
    # make sure the files are not compressed
    with input, util.ungzipped(input.files, tmpdir=tmpdir) as reads:
        seqs, clusterfile = cdhit(input=reads, output=cdhit_tempout,
                                  **cdhit_options)
        # parse raw cd-hit clusters and write it into output_
        with open(clusterfile) as cluster_handle, open(output, 'w') as out:
            for cluster in parse_cdhit_clusters(drop_empty, cluster_handle):
                print('\t'.join(cluster), file=out)
        # delete temporary cd-hit files
        os.remove(seqs)
        os.remove(clusterfile)
        # a specified output destination means that output files can be observed
        # by the callee and their destruction should not be subject to any
        # race conditions
        return data.SampleClusters(input.name, clusters=output,
                                   delete=outdir is None)


@util.fallible(RuntimeError, FileNotFoundError)
def cdpick_multiple(tmpdir: str, input: data.MultipleSampleReads,
                    outdir: Optional[str], drop_empty: bool, **cdhit_options) \
        -> Optional[data.MultipleSampleClusters]:

    return data.MultipleSampleClusters([
        cdpick(tmpdir, sample, outdir, drop_empty, **cdhit_options)
        for sample in input.samples
    ])


if __name__ == '__main__':
    raise RuntimeError