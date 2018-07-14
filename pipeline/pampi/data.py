from typing import Optional, Union, Tuple, Callable, Sequence, Iterable, TypeVar, List, NamedTuple
from contextlib import AbstractContextManager, suppress
from itertools import filterfalse
import abc
import os

from fn import F


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


class SampleReads(SampleFiles):

    def __init__(self, name: str, forward: str, reverse: Optional[str]=None,
                 delete=True):
        reads = (forward, reverse) if reverse else (forward,)
        super().__init__(name, *reads, delete=delete)
        self._paired = reverse is not None

    @property
    def paired(self) -> bool:
        return self._paired


class SampleClusters(SampleFiles):

    def __init__(self, name: str, clusters: str, delete=True):
        super().__init__(name, clusters, delete=delete)

    @property
    def clusters(self) -> Optional[str]:
        return self.files[0] if self.files else None


# TODO we might want to implement full-blown classes with init-time validation
# to make sure MultipleSample* can't be initialised with released resources
MultipleSampleReads = NamedTuple('MultipleSampleReads', [
    ('samples', List[SampleReads])
])

MultipleSampleClusters = NamedTuple('MultipleSampleClusters', [
    ('samples', List[Optional[SampleClusters]])
])

if __name__ == '__main__':
    raise RuntimeError
