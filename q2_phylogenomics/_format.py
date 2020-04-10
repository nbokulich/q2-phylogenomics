# ----------------------------------------------------------------------------
# Copyright (c) 2020, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------


import qiime2.plugin.model as model

import subprocess

from q2_types.feature_data import DNAFASTAFormat


# http://www.htslib.org/doc/samtools-mpileup.html#Pileup_Format
class PileUpTSVFormat(model.TextFileFormat):
    def _validate_(self, level):
        pass


class PileUpFilesDirFmt(model.DirectoryFormat):
    pileups = model.FileCollection(r'.+\.tsv', format=PileUpTSVFormat)

    @pileups.set_path_maker
    def pileups_path_maker(self, name):
        return name + '.tsv'


class GenBankFormat(model.TextFileFormat):
    def _validate_(self, level):
        pass


GenBankDirFmt = model.SingleFileDirectoryFormat(
    'GenBankDirFmt', 'sequence.gb', GenBankFormat)


class BAMFormat(model.BinaryFileFormat):
    def _validate_(self, level):
        cmd = ['samtools', 'quickcheck', '-v', str(self)]
        result = subprocess.run(cmd)
        if result.returncode != 0:
            raise model.ValidationError(
                'samtools quickcheck -v failed on %s' % self.path.name)


class SAMFormat(model.TextFileFormat):
    def _validate_(self, level):
        cmd = ['samtools', 'quickcheck', '-v', str(self)]
        result = subprocess.run(cmd)
        if result.returncode != 0:
            raise model.ValidationError(
                'samtools quickcheck -v failed on %s' % self.path.name)


class BAMFilesDirFmt(model.DirectoryFormat):
    bams = model.FileCollection(r'.+\.bam', format=BAMFormat)

    @bams.set_path_maker
    def bams_path_maker(self, name):
        return name + '.bam'


class SAMFilesDirFmt(model.DirectoryFormat):
    sams = model.FileCollection(r'.+\.sam', format=SAMFormat)

    @sams.set_path_maker
    def sams_path_maker(self, name):
        return name + '.sam'


class FASTAFilesDirFmt(model.DirectoryFormat):
    fastas = model.FileCollection(r'.+\.fasta', format=DNAFASTAFormat)

    @fastas.set_path_maker
    def fastas_path_maker(self, name):
        return name + '.fasta'
