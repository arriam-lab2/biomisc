import os
from typing import Iterable, Tuple, Optional

import numba as nb
import numpy as np
from fn import F

from pipeline.pampi import util, data


@nb.jit(locals={'total': nb.int32, 'threshold': nb.int32, 'stop': nb.int32})
def qualstop(minqual: int, window: int, scores: np.ndarray):
    """
    Return the right non-inclusive border passing the quality threshold
    :param minqual:
    :param window:
    :param scores:
    :return:
    """
    if minqual < 1 or window < 1:
        raise ValueError('minqual and window must be positive')
    if len(scores) < window:
        raise ValueError('window > len(scores)')
    total = scores[:window].sum()
    threshold = minqual * window
    if total < threshold:
        return 0
    stop = window
    for i in range(window, len(scores)):
        total += (scores[i] - scores[i-window])
        if total < threshold:
            break
        stop += 1
    return stop


def decode(base: int, reads: Iterable[Tuple[str, str, str]]) \
        -> Iterable[Tuple[str, str, np.ndarray]]:
    """
    Decode quality strings in raw fastq reads
    :param base: Phred base quality
    :param reads:
    :return:
    """
    # TODO should we perform overflow checks?
    for name, seq, qual in reads:
        decoded = np.frombuffer(qual.encode(), dtype=np.uint8).astype(np.int32)
        yield name, seq, decoded - base


def encode(base: int, reads: Iterable[Tuple[str, str, np.ndarray]]) \
        -> Iterable[Tuple[str, str, str]]:
    """
    Encode quality score arrays in fastq reads
    :param base: Phred base quality
    :param reads:
    :return:
    """
    # TODO should we perform overflow checks?
    for name, seq, qual in reads:
        encoded = (qual+base).astype(np.uint8).tobytes().decode('utf8')
        yield name, seq, encoded


def rolling(minqual, window, reads: Iterable[Tuple[str, str, np.ndarray]]) \
        -> Iterable[Tuple[str, str, np.ndarray]]:
    for name, seq, qual in reads:
        stop = qualstop(minqual, window, qual)
        yield name, seq[:stop], qual[:stop].copy()


def headcrop(length, reads: Iterable[Tuple[str, str, np.ndarray]]) \
        -> Iterable[Tuple[str, str, np.ndarray]]:
    for name, seq, qual in reads:
        yield name, seq[length:], qual[length:].copy()


def trim(phred: int, minqual: int, window: int, minlen: int, croplen: int,
         reads: Iterable[Tuple[str, str, str]]) -> Iterable[Tuple[str, str, str]]:
    return (
        F(decode, phred) >>
        (rolling, minqual, window) >>
        (filter, lambda x: len(x[1]) >= minlen) >>
        (headcrop, croplen) >>
        (encode, phred)
    )(reads)


def cumlength(pair: Tuple[Tuple[str, str, str], Tuple[str, str, str]]) -> int:
    return len(pair[0][1]) + len(pair[1][1])


def trimmer(tmpdir: str, phred: int, minqual: int, window: int, minlen: int,
            croplen: int, compress: bool, outdir: Optional[str],
            samples: data.MultiplePairedFastq) -> data.MultiplePairedFastq:
    # do not filter individual reads by length
    # TODO devectorise all low-level generators to make trimmer_ atomic: ...
    # TODO ... this might shave off redundant complexity below
    trimmer_ = F(trim, phred, minqual, window, 0, croplen)
    trimmed_samples = []
    fwd_suffix = f'_R1.{util.FASTQ}' + util.ending(compress)
    rev_suffix = f'_R2.{util.FASTQ}' + util.ending(compress)

    for sample in samples.samples:
        fwd_out = (
            util.randname(tmpdir, fwd_suffix) if outdir is None else
            os.path.join(outdir, sample.name+fwd_suffix)
        )
        rev_out = (
            util.randname(tmpdir, rev_suffix) if outdir is None else
            os.path.join(outdir, sample.name + rev_suffix)
        )
        # filter pairs with insufficient cumulative length
        trimmed_pairs = (
            F(util.starapply, zip) >>
            (map, trimmer_) >>
            (util.starapply, zip) >>
            (filter, lambda pair: cumlength(pair) >= (minlen - croplen))
        )(sample.parse())
        with util.writer(compress, fwd_out) as fbuffer, util.writer(compress, rev_out) as rbuffer:
            for (fname, fseq, fqual), (rname, rseq, rqual) in trimmed_pairs:
                print('@'+fname, fseq, '+', fqual, sep='\n', file=fbuffer)
                print('@'+rname, rseq, '+', rqual, sep='\n', file=rbuffer)
        trimmed_samples.append(
            data.SamplePairedFastq(sample.name, fwd_out, rev_out, outdir is None)
        )
    return data.MultiplePairedFastq(trimmed_samples)

# TODO implement trimmer for multiple single-end FASTQ files and single samples

if __name__ == '__main__':
    raise RuntimeError
