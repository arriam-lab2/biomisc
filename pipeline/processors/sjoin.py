#! /usr/bin/env python

"""

Join FASTA/FASTQ files into a single dataset. Consider sample files S1.fastq
and S2.fastq wherein reads are not named after the samples. If you plan to do
something with both samples at once, e.g. sequence mapping versus a reference
database, whilst retaining read-sample references, you would want to rename all
reads after the samples they come from. This tool provides two ways to do that.
You can either pass a list of files as arguments to this tool along with a
Python-style regex to extract sample names directly from file basenames or you
can provide a table with desired sample names and files. The naming scheme also
keeps a consistent read counter, thus given the index i of the first renamed
read from sample S2, read S2_j maps directly to the (j-i)-th read in S2.fastq.

"""

import sys

from pipeline import util

if sys.version_info < (3, 6):
    print('This tool requires Python >= 3.6')
    sys.exit(1)

import gzip
import operator as op
import os
import re
from itertools import chain
from typing import Iterator, Iterable, Callable, Mapping, Tuple

try:
    import click
    import tqdm
    from Bio.SeqIO.FastaIO import SimpleFastaParser
    from Bio.SeqIO.QualityIO import FastqGeneralIterator
    from fn import F
except ImportError as err:
    print(f'Missing dependencies: {err}')
    sys.exit(1)

FASTA = 'fasta'
FASTQ = 'fastq'


class BadMapping(ValueError):
    pass


class BadSample(ValueError):
    pass


def parse_mapping(path: str) -> Mapping[str, str]:
    try:
        with open(path) as lines:
            pairs = (
                F(map, str.strip) >> (filter, bool) >>
                (map, lambda x: x.split('\t'))
            )(lines)
            return dict(pairs)
    except (TypeError, ValueError):
        raise BadMapping(f'bad sample mapping {path}')


def build_mapping(nameext: Callable[[str], str], files) -> Mapping[str, str]:
    names = map(nameext, files)
    return dict(zip(names, files))


def make_extractor(pattern: str, group: bool) -> Callable[[str], str]:
    ext = re.compile(pattern).findall if group else re.compile(pattern).split
    return F(os.path.basename) >> ext >> op.itemgetter(0)


def fparse(format_, path) -> Iterator[Tuple[str, str]]:
    """
    Parse a fasta/fastq file and return a generator of (name, sequence) pairs.
    The file can be gzipped. Compression is inferred from the file content.
    :param format_:
    :param path:
    :return:
    """
    try:
        parser = FastqGeneralIterator if format_ == FASTQ else SimpleFastaParser
        with (gzip.open(path, 'rt') if util.isgzipped(path) else open(path)) as buffer:
            yield from parser(buffer)
    except (TypeError, ValueError):
        raise BadSample(f'bad {format_} file {path}')
    except FileNotFoundError:
        raise BadSample(f'sample file not found {path}')


def write(handle, format_, records: Iterable[Tuple]):
    if format_ == FASTA:
        for i, (name, seq, *_) in enumerate(records):
            print(f'>{name}_{i+1}', seq, sep='\n', file=handle)
    elif format_ == FASTQ:
        for i, (name, seq, qual) in enumerate(records):
            print(f'@{name}_{i+1}', seq, '+', qual, sep='\n', file=handle)
    else:
        raise ValueError(f'Unsupported format {format_}')


@click.command('sjoin', help=__doc__,
               context_settings=dict(help_option_names=['-h', '--help']))
@click.option('-m', '--mapping',
              type=click.Path(exists=True, dir_okay=False, resolve_path=True),
              help='A tab-seperated headerless table. The first column '
                   'specifies sample names, while the second column provides '
                   'sample file paths. Relative paths are resolved with '
                   'respect to your current PWD, hence it might be safer to '
                   'use absolute paths.')
@click.option('-p', '--pattern', type=str,
              help='A Python regular expression. By default the expression '
                   'is used to split basenames and select the first value in '
                   'the split. You might instead want to use a group regex to '
                   'extract the first occurrence of the group by specifying '
                   'the --group flag')
@click.option('--group', is_flag=True, default=False)
@click.option('-f', '--format', required=True,
              type=click.Choice([FASTQ, FASTA]))
@click.option('-of', '--output_format', type=click.Choice([FASTQ, FASTA]),
              default=FASTA,
              help='Output format. This option defaults to FASTA, because '
                   'FASTQ output is only compatible with FASTQ input.')
@click.option('-o', '--output',
              type=click.Path(exists=False, dir_okay=False, resolve_path=True),
              help='Destination file path; the tool will stream to stdout if '
                   'none provided. If the path ends with .gz the output will '
                   'be gzipped.')
@click.argument('files', nargs=-1,
                type=click.Path(exists=True, dir_okay=False, resolve_path=True))
def sjoin(mapping, pattern, group, format, output_format, output, files):
    if mapping and (pattern or group):
        raise click.BadOptionUsage(
            "You can't specify --mapping and --pattern at the same time."
        )
    if mapping and files:
        raise click.BadOptionUsage(
            "You can't use --mapping and provide samples as arguments at the "
            "same time."
        )
    if pattern and not files:
        raise click.BadArgumentUsage(
            "No sample files provided."
        )
    if not (mapping or pattern):
        raise click.BadOptionUsage(
            "No sample name resolution options provided, i.e. no --pattern or "
            "--mapping."
        )
    if output_format == FASTQ and format != FASTQ:
        raise click.BadOptionUsage(
            "FASTQ output option is only compatible with FASTQ input."
        )
    try:
        samples = (
            parse_mapping(mapping) if mapping else
            build_mapping(make_extractor(pattern, group), files)
        )
    except BadMapping as err:
        raise click.BadParameter(err)
    reads = chain.from_iterable(
        ((name, seq, *qual) for _, seq, *qual in records) for name, records in
        zip(tqdm.tqdm(samples.keys()), map(F(fparse, format), samples.values()))
    )
    try:
        if not output:
            write(sys.stdout, output_format, reads)
        else:
            compress = output.endswith('.gz')
            with gzip.open(output, 'wt') if compress else open(output, 'w') as out:
                write(out, output_format, reads)
    except BadSample as err:
        raise click.BadParameter(err)


if __name__ == '__main__':
    sjoin()
