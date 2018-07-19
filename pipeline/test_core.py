import itertools
import random
import operator as op
from functools import reduce
from hypothesis import given
from hypothesis.strategies import text, lists, sampled_from
from typing import TypeVar, Sequence, Callable, Generator, List, Tuple
from fn import F
from .core import _starapply


A = TypeVar('A')
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


if __name__ == "__main__":
    raise RuntimeError