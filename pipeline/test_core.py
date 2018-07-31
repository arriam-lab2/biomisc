import itertools
import random
from functools import reduce
from hypothesis import given
from hypothesis.strategies import sampled_from
from typing import TypeVar, Sequence, Callable, Generator, List, \
    Tuple, Type, Generic
from fn import F
from core import _starapply, Map, Router


A = TypeVar('A')
B = TypeVar('B')
random.seed(42)


def _sample(size: int, population: Sequence[A]) -> Sequence[A]:
    return [random.choice(population) for _ in range(size)]


def _generate_int(min_value: int, max_value: int) -> Generator[int, None, None]:
    while True:
        yield random.randrange(min_value, max_value)


def _generate_str(alphabet: List[str], size: int) -> Generator[str, None, None]:
    while True:
        yield "".join(random.choice(alphabet) for _ in range(size))


def _transform_int(arg: int) -> int:
    return hash(arg)


def _transform_str(arg: str) -> int:
    return hash(str)


def _func_template(*args) -> int:
    return sum(t(v) for t, v in zip(*args))


def _generate_starapply_test_data(size: int) -> Tuple[Callable, List[Callable]]:
    MIN_INT_VALUE = -128
    MAX_INT_VALUE = 128
    STRING_SIZE = 10
    ALPHABET = list(map(chr, range(ord("A"), ord("z"))))
    typedefs = _sample(size, 
        [(F(_generate_int, MIN_INT_VALUE, MAX_INT_VALUE)(), _transform_int),
         (F(_generate_str, ALPHABET, STRING_SIZE)(), _transform_str)]
    )
    generators, transforms = list(zip(*typedefs))
    return F(_func_template, transforms), generators 


def _generate_starapply_test(count: int) -> List[Tuple[Callable, List[Callable]]]:
    MAX_ARGS = 10
    return [_generate_starapply_test_data(random.randint(1, MAX_ARGS)) \
            for _ in range(count)]


@given(sampled_from(_generate_starapply_test(100)))
def test_starapply(test_data: Tuple[Callable, List[Callable]]):
    function, generators = test_data
    args = [[next(g) for g in generators]]
    assert reduce(F, itertools.chain([
        function], args))() == _starapply(function, args)


class TypeChecker(Generic[A]):
    def __init__(self, value: A) -> None:
        assert type(value) == A, f"{value} is not of type {A.__repr__()}"


class ValueStorage(TypeChecker[A]):
    def __init__(self, value: A) -> None:
        super().__init__(value)
        self.value = value


class FloatConvertible(ValueStorage[A]):
    def __init__(self, value: B) -> None:
        super().__init__(value)

    def __float__(self) -> float:
        return float(self.value)


class StrConvertible(ValueStorage[A]):
    def __init__(self, value: B) -> None:
        super().__init__(value)
    
    def __str__(self) -> str:
        return str(self.value)


class IntConvertible(ValueStorage[A]):
    def __init__(self, value: A) -> None:
        super().__init__(value)

    def __int__(self) -> int:
        return int(self.value)


class Float(FloatConvertible[float],
            StrConvertible[float],
            IntConvertible[float]):
    def __init__(self, value: float) -> None:
        assert type(value)
        super(Float, self).__init__(value)


class Str(FloatConvertible[str],
          StrConvertible[str],
          IntConvertible[str]):
    def __init__(self, value: str) -> None:
        super(Str, self).__init__(value)


class Int(FloatConvertible[int],
          StrConvertible[int],
          IntConvertible[int]):
    def __init__(self, value: int) -> None:
        super(Int, self).__init__(value)


domains = [int, str, float, Int, Str, Float]


def _generate_map(domain: Type[A], codomain: Type[B]) -> Map[A, B]:
    assert domain in domains
    assert codomain in domains
    return Map(domain, codomain, codomain)


def _generate_router() -> Router:
    pass


if __name__ == "__main__":
    #raise RuntimeError

    m = _generate_map(int, float)
    Router(_generate_str(), [])