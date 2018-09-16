import sys
from contextlib import ExitStack
from itertools import starmap
from typing import Pattern, Tuple, List, Optional, TypeVar, Iterator

import click
import regex as re
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from fn import F

from pipeline.util import gzread, gzwrite

ALPHABET = {
    'A': 'A',
    'C': 'C',
    'G': 'G',
    'T': 'T',
    'R': '[AG]',
    'Y': '[CT]',
    'S': '[GC]',
    'W': '[AT]',
    'K': '[GT]',
    'M': '[AC]',
    'B': '[CGT]',
    'D': '[AGT]',
    'H': '[ACT]',
    'V': '[ACG]',
    'N': '[ACGT]'
}

A = TypeVar('A')


class Seq:
    """
    BioPython's high-level SeqRecord, Seq and SeqIO interfaces are too slow.
    This class provides a similar, but thinner, abstraction over sequences
    """
    def __init__(self, name: str, seq: str, qual=None):
        self.name = name
        self.seq = seq
        self.qual = qual

    def __getitem__(self, item: slice):
        seq = self.seq[item]
        qual = self.qual[item] if self.qual else None
        return type(self)(self.name, seq, qual)


def mkprimer(substitutions: int, primer: str):
    try:
        base = ''.join(ALPHABET[base] for base in primer)
        fuzzy = f'{{s<={substitutions}}}' if substitutions else ''
        return re.compile(f'^(:?{base}){fuzzy}', flags=re.BESTMATCH)
    except KeyError as err:
        raise ValueError(f'unknown base: {err}')


def match(primers: List[Tuple[A, Pattern]], seq: Seq) -> Optional[Tuple[A, Seq]]:
    for flag, primer in primers:
        match_ = primer.match(seq.seq)
        if match_:
            return flag, seq[match_.end():]
    return None


def normalise_pairs(forward, reverse, reads1: Iterator, reads2: Iterator) -> Iterator[Tuple[Seq, Seq]]:
    for r1, r2 in zip(reads1, reads2):
        match1 = match([('F', forward), ('R', reverse)], r1)
        match2 = match([('R', reverse), ('F', forward)], r2)
        # both matches must be positive and come from different primers
        if not (match1 and match2) or match1[0] == match2[0] or match1[0] == 'R':
            yield None
            continue
        yield match1[1], match2[1]


@click.command('primercut')
@click.option('-f', '--forward', type=str, required=True,
              help='IUPAC-encoded forward primer')
@click.option('-r', '--reverse', type=str, required=True,
              help='IUPAC-encoded reverse primer')
@click.option('-m', '--mismatches', type=int, default=0)
@click.argument('inputs', nargs=2, required=True,
                type=click.Path(exists=True, dir_okay=False, resolve_path=True))
@click.argument('outputs', nargs=2, required=True,
                type=click.Path(exists=False, dir_okay=False, resolve_path=True))
def primercut(forward, reverse, mismatches, inputs: Tuple[str, str], outputs: Tuple[str, str]):
    forward, reverse = map(F(mkprimer, mismatches), [forward, reverse])
    in1, in2 = inputs
    out1, out2 = outputs
    bad_pairs = 0
    total_pairs = 0
    template = '@{}\n{}\n+\n{}'
    with ExitStack() as context:
        reads1, reads2 = (
            F(map, gzread) >> (map, context.enter_context) >>
            (map, FastqGeneralIterator) >> (map, F(starmap, Seq))
        )([in1, in2])
        output1, output2 = (
            F(map, gzwrite) >> (map, context.enter_context)
        )([out1, out2])
        normalised_pairs = normalise_pairs(forward, reverse, reads1, reads2)
        for entry in normalised_pairs:
            total_pairs += 1
            if not entry:
                bad_pairs += 1
                continue
            r1, r2 = entry
            print(template.format(r1.name, r1.seq, r1.qual), file=output1)
            print(template.format(r2.name, r2.seq, r2.qual), file=output2)
    good_pairs = total_pairs - bad_pairs
    print(f'Successfully normalised {good_pairs} ({good_pairs/total_pairs:.1%})'
          f' pairs out of {total_pairs}', file=sys.stderr)


if __name__ == '__main__':
    primercut()
