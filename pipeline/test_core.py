import itertools
import random
import string
import copy
from functools import reduce
from fn import F
from queue import Queue
from enum import Enum
from collections import defaultdict
from hypothesis import given
from hypothesis.strategies import sampled_from
from abc import abstractmethod
from multipledispatch import dispatch
from typing import TypeVar, Generator, List, Tuple, \
    Type, Generic, SupportsFloat, SupportsInt, Callable, \
    Any, Union, Sequence
from .core import _starapply, Map, Router, pcompile, \
    AmbiguousError, NoRouteError, RedundancyError


A = TypeVar('A')
B = TypeVar('B')
FloatCompatible = Union[SupportsFloat, str, bytes]
IntCompatible = Union[SupportsInt, str, bytes]
TypingProtocol = SupportsFloat.__mro__[1]


def _sample(size: int, 
            population: Sequence[A]) -> Sequence[A]:
    return [random.choice(population) for _ in range(size)]


def _generate_int(min_value: int,
                  max_value: int) -> Generator[int, None, None]:
    while True:
        yield random.randrange(min_value, max_value)


def _generate_str(alphabet: List[str],
                  size: int) -> Generator[str, None, None]:
    while True:
        yield "".join(random.choice(alphabet) for _ in range(size))


def _transform_int(arg: int) -> int:
    min_int_value = -128
    max_int_value = 127
    return arg + random.randint(min_int_value, max_int_value)


def _transform_str(arg: str) -> str:
    string_size = 10
    postfix = ''.join(random.choices(string.ascii_letters + string.digits, 
                                     k=string_size))
    return f"{arg}_{postfix}"


def _transform_float(arg: float) -> float:
    return arg + random.random()


def _func_template(*args: Any) -> int:
    return sum(t(v) for t, v in zip(*args))


class TypedStorage(Generic[A]):
    def __init__(self, value: A) -> None:
        super().__init__()
        self.value: A = value


class FloatCompatibleStorage(TypedStorage[FloatCompatible], SupportsFloat):
    def __init__(self, value: FloatCompatible) -> None:
        super().__init__(value)

    def __float__(self) -> float:
        return float(self.value)


class IntCompatibleStorage(TypedStorage[IntCompatible], SupportsInt):
    def __init__(self, value: IntCompatible) -> None:
        super().__init__(value)

    def __int__(self) -> int:
        return int(self.value)


class SupportsStr(TypingProtocol):
    __slots__ = ()

    @abstractmethod
    def __str__(self) -> str:
        pass


class StrCompatibleStorage(TypedStorage[A], SupportsStr):
    def __init__(self, value: A) -> None:
        super().__init__(value)

    def __str__(self) -> str:
        return str(self.value)


# possible casts:
# Int -> Int | Float | Str
# Float -> Float | Str
# Str -> Str | Foo
# Foo -> Str

class Float(FloatCompatibleStorage, StrCompatibleStorage):
    def __init__(self, value: FloatCompatible) -> None:
        super().__init__(float(value))


class Int(IntCompatibleStorage, FloatCompatibleStorage, StrCompatibleStorage):
    def __init__(self, value: IntCompatible) -> None:
        super().__init__(int(value))


class Str(StrCompatibleStorage):
    def __init__(self, value: StrCompatibleStorage) -> None:
        super().__init__(str(value))


class Foo:
    def __init__(self, value: SupportsStr) -> None:
        self.value: str = str(value)

    def __str__(self) -> str:
        return self.value


@dispatch(str, SupportsStr)
def _map_f(_: str, arg: SupportsStr) -> str:
    return _transform_str(str(arg))


@dispatch(float, SupportsFloat)
def _map_f(_: float, arg: SupportsFloat) -> float:
    return _transform_float(float(arg))


@dispatch(str, SupportsFloat)
def _map_f(_: str, arg: SupportsFloat) -> str:
    return _transform_str(str(float(arg)))


@dispatch(int, SupportsInt)
def _map_f(_: int, arg: SupportsInt) -> int:
    return _transform_int(int(arg))


@dispatch(float, SupportsInt)
def _map_f(_: float, arg: SupportsInt) -> float:
    return _transform_float(float(int(arg)))


@dispatch(str, SupportsInt)
def _map_f(_: str, arg: SupportsInt) -> str:
    return _transform_str(str(int(arg)))


@dispatch(str, Foo)
def _map_f(_: str, arg: Foo) -> str:
    return _transform_str(str(arg.value))


@dispatch(Foo, str)
def _map_f(_: Foo, arg: str) -> Foo:
    return Foo(_transform_str(arg))


@dispatch(int, SupportsStr)
def _casting_map(_: int, arg: SupportsStr) -> int:
    return hash(str(arg))


D_co = TypeVar('D_co', SupportsFloat, SupportsInt, covariant=True)
C_co = TypeVar('C_co', SupportsFloat, SupportsInt, covariant=True)


def _get_map_function(domain: Type[D_co],
                      codomain: Type[C_co]) -> Callable[[D_co], C_co]:
    dispatched = _map_f.dispatch(codomain, domain)
    if not dispatched:
        raise RuntimeError(f"Could not find a mapping "
                           "from {domain} to {codomain}")
    return F(dispatched, None)


def _get_map(domain: Type[D_co],
             codomain: Type[C_co]) -> Map[D_co, C_co]:
    return Map(domain, codomain, _get_map_function(domain, codomain))


class TypeNode(Generic[A]):
    def __init__(self, value: Type[A]) -> None:
        self.value: Type[A] = value
        self.edge_pairs: List[Tuple['TypeNode[B]', Map[A, B]]] = []

    def add_edge(self, destination: 'TypeNode[B]', edge: Map[A, B]) -> None:
        self.edge_pairs.append((destination, edge))

    @property
    def edges(self) -> List[Map[A, B]]:
        return [edge[1] for edge in self.edge_pairs]

    def __repr__(self) -> str:
        return str(self.value)

    def __eq__(self, other: 'TypeNode[B]') -> bool:
        return self.value == other.value

    def __hash__(self) -> int:
        return id(self)


Graph = List[TypeNode[A]]
SortedGraph = List[List[TypeNode]]


def _generate_graph(nodes: List[TypeNode], count: int) -> SortedGraph:
    return [random.sample(copy.deepcopy(nodes),
                          random.randint(1, len(nodes))) for _ in range(count)]


def _generate_routers(domains: List[Type],
                      num_layers: int) -> Tuple[List[Router], SortedGraph]:
    generated = False
    nodes = [TypeNode(domain) for domain in domains]

    routers, layers = [], []
    while not generated:
        layers = _generate_graph(nodes, num_layers)

        for i in range(0, num_layers - 1):
            has_output_edges = False
            edge_pairs = list(itertools.product(layers[i], layers[i + 1]))
            for pair in edge_pairs:
                begin, end = pair
                try:
                    begin.add_edge(end, _get_map(begin.value, end.value))
                    has_output_edges = True
                except RuntimeError:
                    pass

            if not has_output_edges:
                break
        else:
            generated = True

        routers = [Router('',
                   list(reduce(itertools.chain, [n.edges for n in layers[i]])))
                   for i in range(len(layers) - 1)]

    return routers, layers


def _generate_pcompile_test(domains: List[Type],
                            num_layers: int,
                            num_tests) -> List[Tuple[List[Router],
                                                     SortedGraph,
                                                     Type[A],
                                                     Type[B]]]:
    assert num_layers > 0
    routers, layers = _generate_routers(domains, num_layers)
    return [(routers, layers,
             random.choice(routers[0].domains),
             random.choice(routers[-1].codomains))
            for _ in range(num_tests)]


class BfsResult(Enum):
    REACHABLE = 1,
    UNREACHABLE = 2,
    AMBIGUOUS = 3


def _bfs(graph: Graph, input: TypeNode, output: TypeNode) -> BfsResult:
    queue = Queue()
    queue.put(input)
    visited = defaultdict(bool)
    reachable = False
    while not queue.empty():
        node = queue.get()
        if visited[node]:
            return BfsResult.AMBIGUOUS
        elif node is output:
            reachable = True

        visited[node] = True
        for ep in node.edge_pairs:
            queue.put(ep[0])

    return BfsResult.REACHABLE if reachable else BfsResult.UNREACHABLE


DOMAINS = [Int, Float, Str, int, str, float, Foo]


@given(sampled_from(_generate_pcompile_test(DOMAINS, 5, 100)))
def test_pcompile(test_data: Tuple[List[Router],
                  SortedGraph,
                  Type[A],
                  Type[B]]):
    routers, layers, input, output = test_data
    assert (len(routers) > 0) and (len(routers) == len(layers) - 1)
    assert input in routers[0].domains
    assert output in routers[-1].codomains

    input_node = list(filter(lambda n: n.value == input, layers[0]))[0]
    output_node = list(filter(lambda n: n.value == output, layers[-1]))[0]
    result = _bfs(layers, input_node, output_node)
    try:
        pcompile(routers, input, output)
        assert result == BfsResult.REACHABLE
    except AmbiguousError:
        assert result == BfsResult.AMBIGUOUS
    except NoRouteError:
        assert result == BfsResult.UNREACHABLE


def _generate_rshift_test(domains: List[Type],
                          num_tests) -> List[Tuple[List[Router],
                                                   SortedGraph,
                                                   Type[A],
                                                   Type[B]]]:
    return _generate_pcompile_test(domains, 3, num_tests)


@given(sampled_from(_generate_rshift_test(DOMAINS, 50)))
def test_rshift(test_data: Tuple[List[Router],
                                 SortedGraph,
                                 Type[A],
                                 Type[B]]):
    routers, layers, input, output = test_data
    assert(len(routers) == 2)
    left, right = routers
    input_node = list(filter(lambda n: n.value == input, layers[0]))[0]
    output_node = list(filter(lambda n: n.value == output, layers[-1]))[0]
    result = _bfs(layers, input_node, output_node)
    try:
        pcompile([left >> right], input, output)
        assert result == BfsResult.REACHABLE
    except AmbiguousError:
        assert result == BfsResult.AMBIGUOUS
    except NoRouteError:
        assert result == BfsResult.UNREACHABLE
    except RedundancyError:
        pass


def _generate_starapply_test_data(size: int) -> Tuple[Callable, 
                                                      List[Callable]]:
    MIN_INT_VALUE = -128
    MAX_INT_VALUE = 128
    STRING_SIZE = 10
    ALPHABET = list(map(chr, range(ord("A"), ord("z"))))
    typedefs = _sample(size, [
         (F(_generate_int, MIN_INT_VALUE, MAX_INT_VALUE)(), hash),
         (F(_generate_str, ALPHABET, STRING_SIZE)(), hash)
         ]
    )
    generators, transforms = list(zip(*typedefs))
    return F(_func_template, transforms), generators


def _generate_starapply_test(count: int) -> List[Tuple[Callable,
                                                       List[Callable]]]:
    MAX_ARGS = 15
    return [_generate_starapply_test_data(random.randint(1, MAX_ARGS))
            for _ in range(count)]


@given(sampled_from(_generate_starapply_test(1000)))
def test_starapply(test_data: Tuple[Callable, List[Callable]]):
    function, generators = test_data
    args = [[next(g) for g in generators]]
    assert reduce(F, itertools.chain([
        function], args))() == _starapply(function, args)


if __name__ == "__main__":
    raise RuntimeError
