import abc
import os
from contextlib import AbstractContextManager, suppress
from itertools import filterfalse
from typing import Optional, Callable, Sequence, Iterable, TypeVar, List, \
    NamedTuple, Tuple

from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from fn import F

from pipeline.pampi import util

A = TypeVar('A')


_args: Callable[..., Sequence] = lambda *args: args
_true: Callable[[Iterable[Optional[A]]], Iterable[A]] = F(filter, bool)
_missing: Callable[[Iterable[str]], List[str]] = (
    F(_true) >> (filterfalse, os.path.exists) >> list
)
_nonfile: Callable[[Iterable[str]], List[str]] = (
    F(_true) >> (filterfalse, os.path.isfile) >> list
)


class VolatileResource(AbstractContextManager, metaclass=abc.ABCMeta):

    @property
    @abc.abstractmethod
    def released(self) -> bool:
        pass

    @abc.abstractmethod
    def release(self):
        """
        Release the resource
        :return:
        """
        pass

    def __enter__(self) -> 'VolatileResource':
        if self.released:
            raise RuntimeError('accessing a released resource')
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        # TODO do something with the exceptions
        self.release()


# TODO make this class a real Maybe-like monad with an Empty instance
class SampleFiles(VolatileResource):
    """
    The class is Maybe-like
    """

    # TODO we might want to set delete=False or make the argument mandatory
    def __init__(self, name: str, *files, delete=True):
        missing = _missing(files)
        if missing:
            raise FileNotFoundError(f'missing read file(s): {missing}')
        nonfile = _nonfile(files)
        if nonfile:
            raise ValueError(f'the following paths {nonfile} are not files')
        self._name = name
        self._files = files
        self._released = not files
        self._delete = delete

    def __repr__(self):
        return f'{type(self).__name__}({self.name}, ...)'

    def __bool__(self):
        return self.released or not self.files

    @property
    def name(self) -> str:
        return self._name

    @property
    def files(self) -> Optional[Sequence[str]]:
        if self.released:
            return None
        if _missing(self._files):
            raise RuntimeError(f'accessing corrupted {self!r} resource')
        return self._files

    @property
    def released(self) -> bool:
        return self._released

    def release(self):
        if not self.released and self._delete:
            for fname in self.files:
                with suppress(FileNotFoundError):
                    os.remove(fname)
                # TODO maybe we should throw a warning?
        self._released = True


class SamplePairedFastq(SampleFiles):

    def __init__(self, name: str, forward: str, reverse: str,
                 delete=True):
        super().__init__(name, forward, reverse, delete=delete)

    @property
    def forward(self) -> Optional[str]:
        return self.files[0] if self.files else None

    @property
    def reverse(self) -> Optional[str]:
        return self.files[1] if self.files else None

    def parse(self) \
            -> List[Tuple[Tuple[str, str, str], Tuple[str, str, str]]]:
        if self.released:
            raise RuntimeError(f'accessing a released resource {self}')
        with util.gzread(self.forward) as fwd, util.gzread(self.reverse) as rev:
            return list(zip(*map(FastqGeneralIterator, [fwd, rev])))


class SampleFasta(SampleFiles):

    def __init__(self, name: str, sequences: str, delete=True):
        super().__init__(name, sequences, delete=delete)

    @property
    def sequences(self) -> str:
        return self.files[0] if self.files else None

    def parse(self) -> List[Tuple[str, str]]:
        if self.released:
            raise RuntimeError(f'accessing a released resource {self}')
        with util.gzread(self.sequences) as buffer:
            return list(SimpleFastaParser(buffer))


class SampleFastq(SampleFasta):

    def parse(self) -> List[Tuple[str, str, str]]:
        if self.released:
            raise RuntimeError(f'accessing a released resource {self}')
        with util.gzread(self.sequences) as buffer:
            return list(FastqGeneralIterator(buffer))


class SampleClusters(SampleFiles):

    def __init__(self, name: str, clusters: str, delete=True):
        super().__init__(name, clusters, delete=delete)

    @property
    def clusters(self) -> Optional[str]:
        return self.files[0] if self.files else None

    def parse(self) -> List[Tuple[str, List[str]]]:
        if self.released:
            raise RuntimeError(f'accessing a released resource {self}')
        with util.gzread(self.clusters) as buffer:
            return (
                F(map, str.strip) >> (filter, bool) >>
                (map, lambda x: x.split('\t')) >>
                (map, lambda x: (x[0], x[1:])) >> list
            )(buffer)


# TODO we might want to implement full-blown classes with init-time validation
# to make sure MultipleSample* can't be initialised with released resources
MultipleFasta = NamedTuple('MultipleFasta', [
    ('samples', List[SampleFasta])
])

MultipleFastq = NamedTuple('MultipleFastq', [
    ('samples', List[SampleFastq])
])

MultiplePairedFastq = NamedTuple('MultiplePairedFastq', [
    ('samples', List[SamplePairedFastq])
])

MultipleClusters = NamedTuple('MultipleClusters', [
    ('samples', List[Optional[SampleClusters]])
])


if __name__ == '__main__':
    raise RuntimeError
