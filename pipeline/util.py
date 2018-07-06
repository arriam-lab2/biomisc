import os
import subprocess as sp

DEVNULL = open(os.devnull, 'w')
ENV = '/usr/bin/env'
QUITE = dict(stdout=DEVNULL, stderr=sp.STDOUT)


def onpath(executable: str) -> bool:
    """
    Can /usr/bin/env find the executable?
    :param executable: an executable to find
    :return:
    """
    return sp.run([ENV, executable], **QUITE).returncode != 127


def isgzipped(path: str) -> bool:
    with open(path, 'rb') as buffer:
        return buffer.read(3) == b'\x1f\x8b\x08'