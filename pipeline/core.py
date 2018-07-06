import inspect
from functools import reduce
from typing import Callable, TypeVar, Generic, Type, List, Tuple

# declare types

A = TypeVar('A')
B = TypeVar('B')
C = TypeVar('C')


Map = Tuple[Type[A], Type[B], Callable[[A], B]]


class Node:
    def __init__(self, name: str, maps: List[Map]):
        pass

    def domain(self, dtype: Type[A]) -> List[Map[[A], B]]:
        pass

    def codomain(self, dtype: Type[A]) -> List[Map[[B], A]]:
        pass


def compile(nodes: List[Node]) -> Callable:
    pass


if __name__ == '__main__':
    raise RuntimeError
