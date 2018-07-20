import os
import tempfile
from typing import List, Callable, TypeVar

import click
import pandas as pd
from fn import F, _ as X
from fn.func import identity

from pipeline import core
from pipeline.pampi import data, pick, util, join

CLUSTERS = 'clusters'
TMPDIR = 'tmpdir'
FASTQ = 'fastq'
FASTA = 'fasta'
PAIRED_FASTQ = 'paired_fastq'

A = TypeVar('A')
B = TypeVar('B')


_INPUT_DTYPE_DISPATCH = {
    CLUSTERS: (data.SampleClusters, data.MultipleClusters),
    FASTA: (data.SampleFasta, data.MultipleFasta),
    FASTQ: (data.SampleFastq, data.MultipleFastq),
    PAIRED_FASTQ: (data.SamplePairedFastq, data.MultiplePairedFastq)
}


def validate(f: Callable[[A], bool], transform: Callable[[A], B], message: str,
             ctx, param: str, value: A) -> B:
    try:
        transformed = transform(value)
        if not f(value):
            raise click.BadParameter(message, ctx=ctx, param=param)
        return transformed
    except Exception as err:
        raise click.BadParameter(
            f'validator failed with: {err}', ctx=ctx, param=param
        )


_parse_input: Callable[[str], pd.DataFrame] = (
    lambda x: pd.read_csv(x, sep='\t', header=None, dtype=str)
)
# TODO !!!can't specify pe and se libraries in the same file!!!
_input_paths_exist: Callable[[pd.DataFrame], bool] = (
    lambda df: len(df) and df.iloc[:, 1:].applymap(os.path.exists).all().all()
)

# TODO improve documentation of all features
# TODO separate transformation and validation: it's impossible to handle exceptions...
# TODO ... e.g. during input parsing (when a wrong format is specified)

@click.group(chain=True, invoke_without_command=True,
             context_settings=dict(help_option_names=['-h', '--help']))
@click.option('-i', '--input', required=True,
              type=click.Path(exists=True, dir_okay=False, resolve_path=True),
              callback=F(validate,
                         F(_parse_input) >> _input_paths_exist,
                         _parse_input,
                         'not all paths specified in the input mapping exist '
                         'or the mapping is empty'))
@click.option('-d', '--dtype', required=True,
              type=click.Choice([CLUSTERS, FASTA, FASTQ, PAIRED_FASTQ]),
              help='Initial data type')
@click.option('-t', '--tempdir', default=tempfile.gettempdir(),
              type=click.Path(exists=False, dir_okay=True, resolve_path=True),
              callback=F(validate, os.path.isdir, identity,
                         'tempdir is not a directory or does not exist'),
              help='Temporary directory location')
@click.pass_context
def pampi(ctx, input: pd.DataFrame, dtype: str, tempdir: str):
    ctx.obj[TMPDIR] = tempdir


@pampi.resultcallback()
@click.pass_context
def pipeline(ctx, routers: List[core.Router], input: pd.DataFrame, dtype, *_, **__):
    if not routers:
        exit()
    # TODO streamline input conversion
    # convert parsed data into an appropriate data type
    single_t, multiple_t = _INPUT_DTYPE_DISPATCH[dtype]
    try:
        samples = multiple_t([
            single_t(name, *files, delete=False)
            for name, *files in input.itertuples(False)
        ])
    except (TypeError, IndexError):
        raise ValueError(f'input data are not compatible with data type {dtype}')
    output = core.pcompile(routers, multiple_t, None)(samples)


# TODO add validators

@pampi.command('qc')
@click.pass_context
@click.option('-l', '--minlen', type=int)
@click.option('-q', '--minqual', type=str, help='e.g. 18:3')  # TODO document
@click.option('-c', '--head_crop', type=int)
@click.option('-o', '--output', required=True,
              type=click.Path(exists=False, resolve_path=True),
              callback=F(validate, lambda v: not os.path.exists(v), identity,
                         'output destination exists'),
              help='Output destination.')
@click.pass_context
def qc(ctx):
    pass


@pampi.command('join')
@click.pass_context
@click.option('-p', '--pattern', type=str,
              help='A Python regular expression. By default the expression '
                   'is used to split basenames and select the first value in '
                   'the split. You might instead want to use a group regex to '
                   'extract the first occurrence of the group by specifying '
                   'the --group flag')
@click.option('--group', is_flag=True, default=False)
@click.option('--compress', is_flag=True, default=False)
@click.option('-o', '--output',
              type=click.Path(exists=False, dir_okay=True, resolve_path=True),
              callback=F(validate,
                         lambda x: not x or util.root_exists(x),
                         identity,
                         'destination root does not exist'),
              help='Output destination. This should be a regular path for '
                   'single-file output types or a pattern for paired-end '
                   'fastq output. Pattern example: /path/to/output%.fastq - '
                   'here % will be replaced by R1 and R2 for forward and '
                   'reverse reads respectively')
def joiner(ctx, pattern, group, compress, output):
    rename = join.make_extractor(pattern, group) if pattern else identity
    output = output.replace('%', '{}') if output else None
    options = (ctx.obj[TMPDIR], rename, compress, output)
    return core.Router('joiner', [
        core.Map(data.MultipleFasta, data.SampleFasta,
                 lambda samples: join.join(*options, samples)),
        core.Map(data.MultipleFastq, data.SampleFastq,
                 lambda samples: join.join(*options, samples)),
        core.Map(data.MultiplePairedFastq, data.SamplePairedFastq,
                 lambda samples: join.join(*options, samples)),
        core.Map(data.MultipleClusters, data.SampleClusters,
                 lambda samples: join.join(*options, samples))
    ])


@pampi.command('pick')
@click.option('-r', '--reference', required=True,
              type=click.Path(exists=True, dir_okay=False, resolve_path=True),
              help='Reference dataset (FASTA or FASTQ)')
@click.option('-a', '--accurate', is_flag=True, default=False,
              help='Run in accurate mode. By default, a sequence is clustered '
                   'to the first cluster that meet the threshold (fast '
                   'cluster). In this mode the program will cluster it into '
                   'the most similar cluster that meet the threshold. This '
                   'might take several times more runtime, though.')
@click.option('-s', '--similarity', type=float, default=0.97,
              callback=F(validate, lambda v: 0.5 <= v <= 1, identity,
                         'not in [0.5, 1]'),
              help='Sequence similarity cutoff value; a floating point number '
                   'within [0.5, 1].')
@click.option('-t', '--threads', type=int, default=1,
              callback=F(validate, X > 0, identity, 'must be positive'),
              help='The number of CPU threads to use.')
@click.option('-m', '--memory', type=int, default=1000,
              callback=F(validate, X >= 100, identity,
                         'should be at least 100MB'),
              help='Maximum amount of RAM available to CD-HIT (must be at '
                   'least 100MB).')
@click.option('-e', '--drop_empty', is_flag=True, default=False,
              help='delete empty output')
@click.option('-o', '--outdir',
              type=click.Path(exists=True, dir_okay=True, resolve_path=True),
              callback=F(validate, lambda x: not x or os.path.isdir(x), identity,
                         'destination does not exist or is not a directory'),
              help='Output destination.')
@click.pass_context
def picker(ctx, reference: str, accurate: bool, similarity: float, threads: int,
           memory: int, drop_empty: bool, outdir: str):
    options = dict(tmpdir=ctx.obj[TMPDIR], outdir=outdir, drop_empty=drop_empty,
                   reference=reference, accurate=accurate,
                   similarity=similarity, threads=threads, memory=memory)
    return core.Router('picker', [
        core.Map(data.SampleFasta, data.SampleClusters,
                 lambda x: pick.cdpick(sample=x, **options)),
        core.Map(data.MultipleFasta, data.MultipleClusters,
                 lambda x: pick.cdpick_multiple(samples=x, **options)),
        core.Map(data.SampleFastq, data.SampleClusters,
                 lambda x: pick.cdpick(sample=x, **options)),
        core.Map(data.MultipleFastq, data.MultipleClusters,
                 lambda x: pick.cdpick_multiple(samples=x, **options)),
        core.Map(data.SamplePairedFastq, data.SampleClusters,
                 lambda x: pick.cdpick(sample=x, **options)),
        core.Map(data.MultiplePairedFastq, data.MultipleClusters,
                 lambda x: pick.cdpick_multiple(samples=x, **options))

    ])


if __name__ == '__main__':
    pampi(obj={})
