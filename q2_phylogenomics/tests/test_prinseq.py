# ----------------------------------------------------------------------------
# Copyright (c) 2020, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import gzip
import itertools
import unittest

from q2_types.per_sample_sequences import (
    SingleLanePerSampleSingleEndFastqDirFmt,
    SingleLanePerSamplePairedEndFastqDirFmt,
    FastqGzFormat,
)
from qiime2 import Artifact
from qiime2.plugin.testing import TestPluginBase


class TestPrinseqSingle(TestPluginBase):
    package = 'q2_phylogenomics.tests'

    # This test is really just to make sure that the command runs - the
    # detailed tests in the Util Tests below ensure the commands are crafted
    # appropriately.
    def test_typical(self):
        demuxed_art = Artifact.load(self.get_data_path('single-end.qza'))
        obs_art, = self.plugin.methods['prinseq_single'](demuxed_art)
        obs = obs_art.view(SingleLanePerSampleSingleEndFastqDirFmt)
        obs_seqs = obs.sequences.iter_views(FastqGzFormat)
        for _, obs_fp in obs_seqs:
            with gzip.open(str(obs_fp), 'rt') as obs_fh:
                # Iterate over expected and observed reads, side-by-side
                for records in itertools.zip_longest(*[obs_fh] * 4):
                    (obs_seq_h, obs_seq, _, obs_qual) = records
                    # Make sure low-qual bases were removed
                    self.assertTrue('#' not in obs_qual)
                    # Make sure prinseq trimmed the sequences, too
                    self.assertTrue(len(obs_seq) == len(obs_qual))


class TestPrinseqPaired(TestPluginBase):
    package = 'q2_phylogenomics.tests'

    # This test is really just to make sure that the command runs - the
    # detailed tests in the Util Tests below ensure the commands are crafted
    # appropriately.
    def test_typical(self):
        demuxed_art = Artifact.load(self.get_data_path('paired-end.qza'))
        # The forward and reverse reads are identical in these data
        obs_art, = self.plugin.methods['prinseq_paired'](demuxed_art)
        obs = obs_art.view(SingleLanePerSamplePairedEndFastqDirFmt)
        obs_seqs = obs.sequences.iter_views(FastqGzFormat)
        for _, obs_fp in obs_seqs:
            with gzip.open(str(obs_fp), 'rt') as obs_fh:
                # Iterate over expected and observed reads, side-by-side
                for records in itertools.zip_longest(*[obs_fh] * 4):
                    (obs_seq_h, obs_seq, _, obs_qual) = records
                    # Make sure low-qual bases were removed
                    self.assertTrue('#' not in obs_qual)
                    # Make sure prinseq trimmed the sequences, too
                    self.assertTrue(len(obs_seq) == len(obs_qual))


if __name__ == '__main__':
    unittest.main()
