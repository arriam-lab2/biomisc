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


# external prototype for Processor's compatible
def composable(a: Callable[[A], B], b: Callable[[B], C]) -> bool:
    """
    Make sure functions `a` and `b` can be composed, i.e. (a . b)(x) := b(a(x))
    type-checks. In other words, make sure the functions actually conform to
    this function's type-signature.
    :raises TypeError: if signatures are broken
    >>> def f(x: int) -> str: pass
    >>> def g(x: int) -> None: pass
    >>> def h(x: str) -> int: pass
    >>> composable(f, g)
    False
    >>> composable(f, h)
    True
    >>> composable(g, h)
    False
    >>> composable(h, g)
    True
    """
    rettype = inspect.signature(a).return_annotation
    first_arg, *other = inspect.signature(b).parameters.values()
    argtype = first_arg.annotation
    # output type should be compatible with input type and there should only
    # be one argument
    return issubclass(rettype, argtype) and not other


if __name__ == '__main__':
    raise RuntimeError
