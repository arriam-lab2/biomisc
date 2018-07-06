from typing import List
import operator as op
import tempfile
from functools import reduce

import click

from pipeline import processors, data


# TODO load dtype converters
# TODO load dtype specs dynamically
# TODO load processor specs dynamically


TMPDIR = 'tmpdir'


def validate_input(ctx, param, value) -> data.Samples:
    pass


@click.group(chain=True, invoke_without_command=True)
@click.option('-i', '--input', callback=validate_input, required=True)
@click.option('-t', '--tempdir', default=tempfile.gettempdir())
@click.pass_context
def pampi(ctx, input: data.Samples, tempdir: str):
    ctx.obj[TMPDIR] = tempfile.TemporaryDirectory(dir=tempdir)


@pampi.resultcallback()
@click.pass_context
def pipeline(ctx, processors: List[processors.Processor], input: data.Samples, **_):
    # TODO handle compilation type error
    with ctx.obj[TMPDIR]:
        reduce(op.rshift, processors)(input)


@pampi.command('qc')
@click.pass_context
def qc(ctx):
    prefix: str = ctx.obj[TMPDIR].name
    pass


@pampi.command('join')
@click.pass_context
def join(ctx):
    prefix: str = ctx.obj[TMPDIR].name
    pass


@pampi.command('pick')
@click.pass_context
def pick(ctx):
    prefix: str = ctx.obj[TMPDIR].name
    pass


if __name__ == '__main__':
    pampi(obj={})
