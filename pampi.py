from typing import List
import operator as op
from functools import reduce

import click

from pipeline import processors, data


@click.group(chain=True, invoke_without_command=True)
def cli():
    pass


@cli.resultcallback()
def pipeline(processors: List[processors.Processor], input: data.Samples):
    reduce(op.rshift, processors)(input)


@cli.command('qc')
@click.option('-o', '--option')
def qc(option):
    # TODO load processor
    pass


@cli.command('join')
def join():
    # TODO load processor
    pass


@cli.command('pick')
@click.option('-o', '--option')
def pick(option):
    # load processor
    pass


if __name__ == '__main__':
    cli()
