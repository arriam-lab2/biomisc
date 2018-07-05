from typing import NamedTuple, List, Optional, Union


Reads = NamedTuple('Reads', [
    ('names', List[str]), ('r1', List[str]), ('r2', Optional[List[str]])
])

Mappings = NamedTuple('Tables', [
    ('names', List[str]), ('mappings', List[str])
])


Samples = Union[Reads, Mappings]


if __name__ == '__main__':
    raise RuntimeError
