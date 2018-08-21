import inspect
import operator as op
from functools import reduce
from itertools import combinations, starmap, product, chain
from typing import Callable, TypeVar, Generic, Type, List, Tuple, Optional, \
    Iterable, Sequence, Union

from fn import F

# declare types

A = TypeVar('A')
B = TypeVar('B')
C = TypeVar('C')

POS_ARGS = frozenset(
    [inspect.Parameter.POSITIONAL_ONLY,
     inspect.Parameter.POSITIONAL_OR_KEYWORD,
     inspect.Parameter.VAR_POSITIONAL]
)

_starapply: Callable[[Callable[..., A], Iterable], A] = lambda f, x: f(*x)

# TODO research a way to relax the algebraic type system within Python
# TODO https://github.com/HypothesisWorks/hypothesis/pull/643 (as an example)
# abstract away type compatibility check to not depend on its implementation
_typematch: Union[Callable[[Type, Type], bool]] = op.is_

# signature matches
_signmatch: Callable[[Tuple[Type, Tuple], Tuple[Type, Tuple]], bool] = (
    F(_starapply, zip) >> (starmap, _typematch) >> all
)
# left-side codomain matches right-side domain
_composable: Callable[['Map', 'Map'], bool] = (
    lambda l, r: _typematch(l.codomain, r.domain)
)
# are there any maps with identical domain-codomain pairs?
_redundant: Callable[[Iterable['Map']], bool] = (
        F(map, lambda m: m.signature) >> F(combinations, r=2) >>
        (map, _signmatch) >> any
)


class RedundancyError(ValueError):
    pass


class AmbiguousError(ValueError):
    def __init__(self, input: str, output: str) -> None:
        super().__init__(f'there are several valid routes '
                         'from {input} to {output}')


class NoRouteError(ValueError):
    def __init__(self, input: str, output: str) -> None:
        super().__init__(f'there is no route from {input} to {output}')


class Map(Generic[A, B]):

    def __init__(self, domain: Type[A], codomain: Type[B], f: Callable[[A], B]):
        """
        :param domain: `None` is treated as NoneType
        :param codomain: `None` is treated as NoneType
        :param f:
        """
        # validate domain and codomain
        self._domain = type(None) if domain is None else domain
        self._codomain = type(None) if codomain is None else codomain
        if not inspect.isclass(self._domain):
            raise ValueError(f'domain is not a type')
        if not inspect.isclass(self._codomain):
            raise ValueError(f'codomain is not a type')

        # make sure the function can be called with a single argument
        if not callable(f):
            raise ValueError('f is not callable')
        signature = inspect.signature(
            f if inspect.isfunction(f) else f.__call__
        )
        # leave only positional arguments without default values and *args
        positional = [param for param in signature.parameters.values()
                      if param.kind in POS_ARGS]
        nrequired = sum(p.default is inspect.Parameter.empty for p in positional
                        if p.kind is not inspect.Parameter.VAR_POSITIONAL)
        if nrequired > 1:
            raise ValueError('f has more than one required positional argument')
        if len(positional) < 1:
            raise ValueError('f has no positional arguments')

        self._f: F = f if isinstance(f, F) else F(f)

    @property
    def domain(self) -> Type[A]:
        return self._domain

    @property
    def codomain(self) -> Type[B]:
        return self._codomain

    @property
    def signature(self) -> Tuple[Type[A], Type[B]]:
        return self.domain, self.codomain

    def __repr__(self):
        # TODO maybe we should show show type reprs instead of their names?
        try:
            return f'({self.domain.__name__}) -> {self.codomain.__name__}'
        except AttributeError:
            return f'({self.domain}) -> {self.codomain}'

    def __call__(self, value: A) -> B:
        return self._f(value)

    def __rshift__(self, other: 'Map[B, C]') -> 'Map[A, C]':
        if not isinstance(other, type(self)):
            # TODO maybe we should show type(self) instead of its name?
            raise ValueError(
                f'right-hand operand is not an instance of '
                f'{type(self).__name__}'
            )
        if not _composable(self, other):
            raise ValueError(
                f'domain of right-hand operand {other} does not match codomain '
                f'of {self}'
            )
        return type(self)(self.domain, other.codomain, self._f >> other._f)


class Router:
    """
    A set of computation graph edges.
    """

    def __init__(self, name: str, maps: Sequence[Map]):
        self._name = name
        self._maps = tuple(maps)
        if not all(isinstance(m, Map) for m in maps):
            raise ValueError(f'not all maps are {Map.__name__} instances')
        if _redundant(self._maps):
            raise RedundancyError('mappings are redundant')

    @property
    def name(self) -> str:
        return self._name

    @property
    def maps(self) -> Sequence[Map]:
        return self._maps

    @property
    def domains(self) -> List[Type]:
        """
        Return all available domains
        """
        return [m.domain for m in self._maps]

    @property
    def codomains(self) -> List[Type]:
        """
        Return all available codomains
        """
        return [m.codomain for m in self._maps]

    @property
    def signatures(self) -> List[Tuple[Type, Type]]:
        return [m.signature for m in self._maps]

    def __len__(self):
        return len(self._maps)

    def __bool__(self):
        return bool(len(self))

    def __rshift__(self, other: 'Router') -> 'Router':
        if not isinstance(other, type(self)):
            raise ValueError(
                f'right-hand operand is not an instance of '
                f'{type(self).__name__}'
            )
        name = f'{self.name} -> {other.name}'
        if not (self and other):
            return Router(name, [])
        left = self.constrain(None, other.domains)
        right = other.constrain(left.codomains, None)
        compositions = (
            F(product) >> (filter, F(_starapply, _composable)) >>
            (starmap, op.rshift) >> list
        )(left.maps, right.maps)
        return type(self)(name, compositions)

    def constrain(self,
                  domains: Optional[List[Type]],
                  codomains: Optional[List[Type]]) -> 'Router':
        """
        :param domains: None is regarded as no constraints
        :param codomains: None is regarded as no constraints
        :return:
        """
        if not all(map(inspect.isclass, chain(domains or [], codomains or []))):
            raise ValueError('not all type constraints are types')

        def match_any(t: Type, options: Optional[List[Type]]) -> bool:
            return options is None or any(map(F(_typematch, t), options))

        if domains is None and codomains is None:
            return self

        mappings = [
            m for m in self._maps if
            match_any(m.domain, domains) and match_any(m.codomain, codomains)
        ]
        return type(self)(self.name, mappings)


def pcompile(routers: List[Router], input: Type[A], output: Type[B]) \
        -> Callable[[A], B]:
    """
    Compile a path from input node to the output node
    :param routers:
    :param input:
    :param output:
    :return:
    """
    # TODO maybe we should create and use something like CompileTimeError?
    if not all(isinstance(router, Router) for router in routers):
        raise ValueError(f'not all routers are instances of {Router.__name__}')
    if len(routers) == 1:
        constrained = routers[0].constrain([input], [output])
    else:
        head, *tail = routers
        try:
            composed: Router = reduce(op.rshift, tail, head.constrain([input], None))
        except RedundancyError:
            raise AmbiguousError(input.__name__, output.__name__)
        constrained = composed.constrain(None, [output])
    if not constrained:
        raise NoRouteError(input.__name__, output.__name__)

    return constrained.maps[0]


if __name__ == '__main__':
    raise RuntimeError
