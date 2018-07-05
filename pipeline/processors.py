import inspect
from functools import reduce
from typing import Callable, TypeVar, Generic, Type

# declare types

FASTA = 'FASTA'
FASTQ = 'FASTQ'

A = TypeVar('A')
B = TypeVar('B')
C = TypeVar('C')


# TODO use a proper Maybe monad implementation instead of Optional
# TODO add default implementation to Processor's `compatible` and  `__rshift__`

class Processor(Generic[A, B]):

    def __init__(self, f: Callable[[A], B], argtype: Type, rettype: Type,
                 argtype_validator: Callable[[Type], bool]=None):
        # TODO verify basic signature properties of `f` and `validator`
        # TODO maybe we should only accept callables with proper annotations
        # to avoid explicit argtype/rettype parameters?
        self._functions = (f,)
        self._argtype = argtype
        self._rettype = rettype
        self._validator = (
                argtype_validator or (lambda v: issubclass(v, self._argtype))
        )

    def __str__(self):
        return f'({self._argtype}) -> {self._rettype}'

    def __call__(self, data: A) -> B:
        return reduce(lambda val, f: f(val), self._functions, data)

    def __rshift__(self, acceptor: 'Processor[B, C]') -> 'Processor[A, C]':
        if not acceptor._validator(self._rettype):
            raise TypeError(f'{self} is not compatible with {acceptor}')
        composed = self.__new__(type(self))
        composed._functions = (*self._functions, *acceptor._functions)
        composed._argtype = self._argtype
        composed._rettype = acceptor._rettype
        composed._validator = self._validator
        return composed


# external prototype for Processor's compatible
def _composable(a: Callable[[A], B], b: Callable[[B], C]) -> bool:
    """
    Make sure functions `a` and `b` can be composed, i.e. (a . b)(x) := b(a(x))
    type-checks. In other words, make sure the functions actually conform to
    this function's type-signature.
    :raises TypeError: if signatures are broken
    >>> def f(x: int) -> str: pass
    >>> def g(x: int) -> None: pass
    >>> def h(x: str) -> int: pass
    >>> _composable(f, g)
    False
    >>> _composable(f, h)
    True
    >>> _composable(g, h)
    False
    >>> _composable(h, g)
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
