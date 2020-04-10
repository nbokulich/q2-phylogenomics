# ----------------------------------------------------------------------------
# Copyright (c) 2020, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from q2_types.feature_data import DNAFASTAFormat
from q2_types.bowtie2 import Bowtie2IndexDirFmt

from .plugin_setup import plugin
from ._format import BAMFilesDirFmt, SAMFilesDirFmt, SAMFormat, BAMFormat
from ._util import run_command


@plugin.register_transformer
def _0(dirfmt: SAMFilesDirFmt) -> BAMFilesDirFmt:
    result = BAMFilesDirFmt()

    for path, view in dirfmt.sams.iter_views(SAMFormat):
        samtools_view_cmd = ['samtools', 'view',
                             '-O', 'BAM', str(dirfmt.path / path),
                             '-o', str(result.path / path.stem) + '.bam']
        run_command(samtools_view_cmd)

    return result


@plugin.register_transformer
def _1(dirfmt: BAMFilesDirFmt) -> SAMFilesDirFmt:
    result = SAMFilesDirFmt()

    for path, view in dirfmt.bams.iter_views(BAMFormat):
        samtools_view_cmd = ['samtools', 'view', '-h',
                             '-O', 'SAM', str(dirfmt.path / path),
                             '-o', str(result.path / path.stem) + '.sam']
        run_command(samtools_view_cmd)

    return result


@plugin.register_transformer
def _2(db: Bowtie2IndexDirFmt) -> DNAFASTAFormat:
    result = DNAFASTAFormat()
    bowtie2_inspect_cmd = ['bowtie2-inspect',
                           str(db.path / db.get_basename())]
    run_command(bowtie2_inspect_cmd,
                stdout=open(str(result), 'w'))
    return result
