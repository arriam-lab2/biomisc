import subprocess as sp
import contextlib
import tempfile
import logging
import shutil
import uuid
import gzip
import os
from itertools import repeat, filterfalse
from functools import wraps
from typing import Callable, TypeVar, Optional, Sequence, TextIO

from fn import F, _ as X


DEVNULL = open(os.devnull, 'w')
ENV = '/usr/bin/env'
GZIP = 'gzip'
QUITE = dict(stdout=DEVNULL, stderr=sp.STDOUT)
# extensions
GZ = 'gz'
FASTQ = 'fastq'
FASTA = 'fasta'
CLUSTERS = 'clstr'


A = TypeVar('A')
NoneType = type(None)


writer: Callable[[bool, str], TextIO] = lambda compress, path: (
    gzip.open(path, 'wt') if compress else open(path, 'w')
)
ending: Callable[[bool], str] = (
    lambda compress: f'.{util.GZ}' if compress else ''
)


def starapply(f, args):
    return f(*args)


def root_exists(path: str) -> bool:
    """
    Does root directory of a path exist?
    :param path:
    :return:
    """
    return os.path.exists(os.path.dirname(path))


def randname(dirname: str, suffix: str, check: bool=True) -> str:
    """
    Generate a random file name.
    at the call-time
    :param dirname: specify a root directory
    :param suffix: file extension
    :param check: ensures that such a name does not exist in the file system
    :return:
    """
    return (
        F(map, str) >> (map, F(os.path.join, dirname)) >> (map, X+suffix) >>
        (filterfalse, (os.path.exists if check else lambda x: False)) >> next
    )(uuid.uuid4() for _ in repeat(None))


# TODO replace this temporary lifting solution for Optionals with proper monads
def fallible(*exceptions, logger=None) \
        -> Callable[[Callable[..., A]], Callable[..., Optional[A]]]:
    """
    :param exceptions: a list of exceptions to catch
    :param logger: specify a logger to use; None means the default logger,
    False disables logging altogether.
    """
    def fwrap(f: Callable[..., A]) -> Callable[..., Optional[A]]:

        @wraps(f)
        def wrapped(*args, **kwargs):
            try:
                return f(*args, **kwargs)
            except exceptions:
                message = f'called {f} with *args={args} and **kwargs={kwargs}'
                if logger:
                    logger.exception(message)
                if logger is None:
                    logging.exception(message)
                return None

        return wrapped

    return fwrap


def isgzipped(path: str) -> bool:
    with open(path, 'rb') as buffer:
        return buffer.read(3) == b'\x1f\x8b\x08'


def gzread(path: str) -> TextIO:
    return gzip.open(path, 'rt') if isgzipped(path) else open(path)


@contextlib.contextmanager
def ungzipped(*paths, tmpdir=tempfile.gettempdir()) -> Sequence[str]:
    """
    Takes a sequence of paths (p1, p2, ..., pn) and returns (d1, d2, ..., dn),
    where pi == di if not isgzipped(pi), otherwise di point to a temporary
    decompressed file which will be erased upon exit from the context manager.
    :param paths:
    :param tmpdir: a location for temporary decompressed files
    :return: a sequence of file path strings
    """
    # TODO consider defaulting to the built-in gzip in this case
    gzip_exec = shutil.which(GZIP)
    if not gzip_exec:
        raise RuntimeError(
            'No gzip executable found; is it on your PATH?'
        )
    if not os.path.exists(tmpdir):
        raise ValueError(
            f'directory {tmpdir}, specified as tmpdir, does not exist'
        )
    uncompressed = []
    with contextlib.ExitStack() as stack:
        for path in paths:
            # if gzipped, make a temporary decompressed file
            # shelling gzip out, because Python's built-in gzip implementation
            # is notoriously slow
            # !!! note: opening several connections to a named temporary file
            #           is only possible on Unix-like systems. This function is
            #           thus not Windows-friendly.
            if isgzipped(path):
                buffer = tempfile.NamedTemporaryFile(dir=tmpdir)
                stack.enter_context(buffer)
                sp.run([gzip_exec, '-cdf', path], stdout=buffer)
                buffer.flush()
                uncompressed.append(buffer.name)
            else:
                uncompressed.append(path)
        yield tuple(uncompressed)


if __name__ == '__main__':
    raise RuntimeError
