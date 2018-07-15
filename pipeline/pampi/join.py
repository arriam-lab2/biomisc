import gzip
import operator as op
import os
import re
from typing import Callable, Iterator, Tuple, Iterable

from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from fn import F

from pipeline.pampi import util

FASTA = 'fasta'
FASTQ = 'fastq'


class BadSample(ValueError):
    pass


def make_extractor(pattern: str, group: bool) -> Callable[[str], str]:
    ext = re.compile(pattern).findall if group else re.compile(pattern).split
    return F(os.path.basename) >> ext >> op.itemgetter(0)


def parse(format_, path) -> Iterator[Tuple[str, str]]:
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


if __name__ == '__main__':
    raise RuntimeError
