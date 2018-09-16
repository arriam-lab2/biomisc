import operator as op
import os
import re
from contextlib import ExitStack
from functools import reduce
from itertools import count, islice, chain, tee
from typing import Callable, Iterator, Tuple, Iterable, List

from fn import F
from fn.iters import group_by
from multipledispatch import dispatch

from pipeline.pampi import data
from pipeline import util
from pipeline.util import root_exists, writer, ending


class BadSample(ValueError):
    pass


def make_extractor(pattern: str, group: bool) -> Callable[[str], str]:
    # TODO add docs
    extractor = (re.compile(pattern).findall if group else
                 re.compile(pattern).split)
    return F(os.path.basename) >> extractor >> op.itemgetter(0)

# open sample and read clusters
# parse cluster files as mappings
# rename reads in each mapping
# merge reads into one single file
# select output format: clusters, tsv, biom.


def _make_counter(start: int, template: str) -> Iterator[str]:
    # TODO docs
    return map(template.format, count(start))


def join_fastqc(name_templates: Iterable[str],
                fastqs: Iterable[Iterable[Tuple[str, str, str]]]) \
        -> Iterator[Tuple[str, str, str]]:
    """
    :param name_templates:
    :param fastqs:
    :return:
    >>> fastqcs = [
    ...     [('name', 'seq1', 'qual1'), ('name', 'seq2', 'qual2')],
    ...     [('name', 'seq1', 'qual1'), ('name', 'seq2', 'qual2'),
    ...      ('name', 'seq3', 'qual3')]
    ... ]
    >>> name_templates = ['sample1_{}', 'sample2_{}']
    >>> joined_fastqc = [
    ...     ('sample1_1', 'seq1', 'qual1'),
    ...     ('sample1_2', 'seq2', 'qual2'),
    ...     ('sample2_1', 'seq1', 'qual1'),
    ...     ('sample2_2', 'seq2', 'qual2'),
    ...     ('sample2_3', 'seq3', 'qual3')
    ... ]
    >>> list(join_fastqc(name_templates, fastqcs)) == joined_fastqc
    True
    """
    name_generators = [_make_counter(1, template) for template in name_templates]
    for names, reads in zip(name_generators, fastqs):
        yield from ((name, seq, qual) for name, (_, seq, qual) in zip(names, reads))


def join_fasta(name_templates: Iterable[str],
               fastas: Iterable[Iterable[Tuple[str, str]]]) \
        -> Iterable[Tuple[str, str]]:
    """
    :param name_templates:
    :param fastas:
    :return:
    >>> fastas = [
    ...     [('name', 'seq1'), ('name', 'seq2')],
    ...     [('name', 'seq1'), ('name', 'seq2'), ('name', 'seq3')]
    ... ]
    >>> name_templates = ['sample1_{}', 'sample2_{}']
    >>> joined_fastas = [
    ...     ('sample1_1', 'seq1'),
    ...     ('sample1_2', 'seq2'),
    ...     ('sample2_1', 'seq1'),
    ...     ('sample2_2', 'seq2'),
    ...     ('sample2_3', 'seq3')
    ... ]
    >>> list(join_fasta(name_templates, fastas)) == joined_fastas
    True
    """
    name_generators = [_make_counter(1, template) for template in name_templates]
    for names, reads in zip(name_generators, fastas):
        yield from ((name, seq) for name, (_, seq) in zip(names, reads))


def join_clusters(name_templates: Iterable[str],
                  clusters: Iterable[List[Tuple[str, List]]]) \
        -> List[Tuple[str, List[str]]]:
    """
    :param name_templates:
    :param clusters:
    :return:
    >>> clusters = [
    ...     [
    ...         ('cluster1', [1, 2]),
    ...         ('cluster2', [1, 2])
    ...     ],
    ...     [
    ...         ('cluster1', [1]),
    ...         ('cluster3', [1, 2, 3])
    ...     ]
    ... ]
    >>> name_templates = ['sample1_{}', 'sample2_{}']
    >>> clusters_joined = [
    ...     ('cluster1', ['sample1_1', 'sample1_2', 'sample2_1']),
    ...     ('cluster2', ['sample1_3', 'sample1_4']),
    ...     ('cluster3', ['sample2_2', 'sample2_3', 'sample2_4'])
    ... ]
    >>> join_clusters(name_templates, clusters) == clusters_joined
    True
    """
    name_generators = [_make_counter(1, template) for template in name_templates]
    renamed = (
        ((cls, list(islice(names, 0, len(seqs)))) for cls, seqs in clusters_)
        for names, clusters_ in zip(name_generators, clusters)
    )
    grouped = group_by(op.itemgetter(0), chain.from_iterable(renamed))
    # flatten
    return [
        (key, reduce(op.iadd, map(op.itemgetter(1), group)))
        for key, group in grouped.items()
    ]


@dispatch(str, Callable, bool, (str, type(None)), data.MultipleFastq)
def join(tmpdir: str, rename: Callable[[str], str], compress: bool,
         output: str, samples: data.MultipleFastq) \
        -> data.SampleFastq:

    output_ = (
        output if output is not None else
        util.randname(tmpdir, f'.{util.FASTQ}' + ending(compress))
    )

    if not root_exists(output_):
        raise ValueError(f'missing directory for {os.path.dirname(output_)}')

    with ExitStack() as context:
        samples_: List[data.SampleFastq] = list(
            map(context.enter_context, samples.samples)
        )
        name_templates = [f'@{rename(sample.name)}_{{}}' for sample in samples_]
        reads = (sample.parse() for sample in samples_)
        buffer = context.enter_context(writer(compress, output_))
        for name, seq, qual in join_fastqc(name_templates, reads):
            print(name, seq, '+', qual, sep='\n', file=buffer)
        return data.SampleFastq(
            'joined', output_, output is None
        )


@dispatch(str, Callable, bool, (str, type(None)), data.MultiplePairedFastq)
def join(tmpdir: str, rename: Callable[[str], str], compress: bool,
         output_pattern: str, samples: data.MultiplePairedFastq) \
        -> data.SamplePairedFastq:

    out_pattern = (
        output_pattern if output_pattern is not None else
        util.randname(tmpdir, f'_{{}}.{util.FASTQ}' + ending(compress), False)
    )

    if out_pattern and not root_exists(out_pattern):
        raise ValueError(f'missing directory {os.path.dirname(out_pattern)}')

    fwd_output = out_pattern.format('R1')
    rev_output = out_pattern.format('R2')

    with ExitStack() as context:
        samples_: List[data.SamplePairedFastq] = list(
            map(context.enter_context, samples.samples)
        )
        fwd_name_templates = [
            f'@{rename(sample.name)}_{{}}/1' for sample in samples_
        ]
        rev_name_templates = [
            f'@{rename(sample.name)}_{{}}/2' for sample in samples_
        ]
        fwd_stream, rev_stream = tee((sample.parse() for sample in samples_), 2)
        fwd_reads = join_fastqc(fwd_name_templates,
                                F(map, F(map, op.itemgetter(0)))(fwd_stream))
        rev_reads = join_fastqc(rev_name_templates,
                                F(map, F(map, op.itemgetter(1)))(rev_stream))
        forward_buffer = context.enter_context(writer(compress, fwd_output))
        reverse_buffer = context.enter_context(writer(compress, rev_output))
        for (fname, fseq, fqual), (rname, rseq, rqual) in zip(fwd_reads, rev_reads):
            print(fname, fseq, '+', fqual, sep='\n', file=forward_buffer)
            print(rname, rseq, '+', rqual, sep='\n', file=reverse_buffer)
        return data.SamplePairedFastq(
            'joined', fwd_output, rev_output,
            output_pattern is None
        )


@dispatch(str, Callable, bool, (str, type(None)), data.MultipleFasta)
def join(tmpdir: str, rename: Callable[[str], str], compress: bool,
         output: str, samples: data.MultipleFasta) \
        -> data.SampleFasta:

    output_ = (
        output if output is not None else
        util.randname(tmpdir, f'.{util.FASTA}' + ending(compress))
    )

    if not root_exists(output_):
        raise ValueError(f'missing directory for {os.path.dirname(output_)}')

    with ExitStack() as context:
        samples_: List[data.SampleFasta] = list(
            map(context.enter_context, samples.samples)
        )
        name_templates = [f'>{rename(sample.name)}_{{}}' for sample in samples_]
        reads = (sample.parse() for sample in samples_)
        buffer = context.enter_context(writer(compress, output_))
        for name, seq in join_fasta(name_templates, reads):
            print(name, seq, sep='\n', file=buffer)
        return data.SampleFasta(
            'joined', output_, output is None
        )


@dispatch(str, Callable, bool, (str, type(None)), data.MultipleClusters)
def join(tmpdir: str, rename: Callable[[str], str], compress: bool,
         output: str, samples: data.MultipleClusters) \
        -> data.SampleClusters:

    output_ = (
        output if output is not None else
        util.randname(tmpdir, f'.{util.CLUSTERS}' + ending(compress))
    )

    if not root_exists(output_):
        raise ValueError(f'missing directory for {os.path.dirname(output_)}')

    with ExitStack() as context:
        samples_: List[data.SampleClusters] = list(
            map(context.enter_context, samples.samples)
        )
        name_templates = [f'{rename(sample.name)}_{{}}' for sample in samples_]
        clusters = (sample.parse() for sample in samples_)
        buffer = context.enter_context(writer(compress, output_))
        for name, reads in join_clusters(name_templates, clusters):
            print(name, '\t'.join(reads), sep='\t', file=buffer)
        return data.SampleClusters(
            'joined', output_, output is None
        )


if __name__ == '__main__':
    raise RuntimeError
