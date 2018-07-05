from typing import NamedTuple, Callable, List, Union, TypeVar, Generic, Optional, Iterable, Type
import functools
import abc
import inspect

# declare types

FASTA = 'FASTA'
FASTQ = 'FASTQ'

A = TypeVar('A')
B = TypeVar('B')
C = TypeVar('C')


# TODO use a proper Maybe monad implementation instead of Optional
# TODO add default implementation to Processor's `compatible` and  `__rshift__`

class Processor(Generic[A, B], metaclass=abc.ABCMeta):

    @abc.abstractmethod
    def __repr__(self):
        pass

    @abc.abstractmethod
    def __call__(self, data: A) -> B:
        pass

    @abc.abstractmethod
    def __rshift__(self, other: 'Processor[B, C]') -> 'Processor[A, C]':
        pass

    @abc.abstractmethod
    def compatible(self, data: Type[A]) -> bool:
        """
        Can this processor accept value of Type[A]?
        :param data:
        :return:
        """
        pass



if __name__ == '__main__':
    raise RuntimeError
