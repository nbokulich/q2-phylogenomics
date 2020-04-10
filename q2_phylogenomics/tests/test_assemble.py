# ----------------------------------------------------------------------------
# Copyright (c) 2020, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import unittest

import qiime2
from qiime2.plugin.testing import TestPluginBase

from q2_phylogenomics._format import SAMFormat, SAMFilesDirFmt


class TestMapPaired(TestPluginBase):
    package = 'q2_phylogenomics.tests'

    def setUp(self):
        super().setUp()

        self.demuxed_art = qiime2.Artifact.load(
            self.get_data_path('paired-end.qza'))
        self.paired_mapped_unsorted = qiime2.Artifact.load(
            self.get_data_path('paired-end-mapped-unsorted.qza'))
        self.indexed_genome = qiime2.Artifact.load(
            self.get_data_path('sars2-indexed.qza'))
        self.sorted_alignment_maps = qiime2.Artifact.load(
            self.get_data_path('sorted-alignment-maps.qza'))

    def test_map_paired_mapped_only(self):
        obs_art, = self.plugin.methods['map_paired_reads'](
            self.demuxed_art, self.indexed_genome)
        obs = obs_art.view(SAMFilesDirFmt)
        exp = [('sample_a.sam', 10,
                ('SARS2:6:73:567:7631#0', 'SARS2:6:73:233:3421#0'),
                ('SARS2:6:73:356:9806#0', 'SARS2:6:73:356:9806#0')),
               ('sample_b.sam', 10,
                ('SARS2:6:73:941:1973#0', 'SARS2:6:73:552:2457#0'),
                ('SARS2:6:73:356:9806#0', 'SARS2:6:73:356:9806#0')),
               ('sample_c.sam', 10,
                ('SARS2:6:73:231:3321#0', 'SARS2:6:73:552:2457#0'),
                ('SARS2:6:73:356:9806#0', 'SARS2:6:73:356:9806#0'))]
        obs_sams = obs.sams.iter_views(SAMFormat)
        for (obs_fn, obs_sam), (exp_fn, num_records, mapped_ids, unmapped_ids)\
                in zip(obs_sams, exp):
            self.assertEqual(str(obs_fn), exp_fn)
            with open(str(obs_sam)) as sam_f:
                obs_mapped_ids = [line.split('\t')[0]
                                  for line in sam_f
                                  if not line.startswith('@')]
            self.assertEqual(len(obs_mapped_ids), num_records)
            for e in mapped_ids:
                self.assertIn(e, obs_mapped_ids)
            for e in unmapped_ids:
                self.assertNotIn(e, obs_mapped_ids)

    def test_map_paired_not_mapped_only(self):
        obs_art, = self.plugin.methods['map_paired_reads'](
            self.demuxed_art, self.indexed_genome, mapped_only=False)
        obs = obs_art.view(SAMFilesDirFmt)
        exp = [('sample_a.sam', 12,
                ('SARS2:6:73:567:7631#0', 'SARS2:6:73:233:3421#0'),
                ('SARS2:6:73:356:9806#0', 'SARS2:6:73:356:9806#0')),
               ('sample_b.sam', 12,
                ('SARS2:6:73:941:1973#0', 'SARS2:6:73:552:2457#0'),
                ('SARS2:6:73:356:9806#0', 'SARS2:6:73:356:9806#0')),
               ('sample_c.sam', 12,
                ('SARS2:6:73:231:3321#0', 'SARS2:6:73:552:2457#0'),
                ('SARS2:6:73:356:9806#0', 'SARS2:6:73:356:9806#0'))]
        obs_sams = obs.sams.iter_views(SAMFormat)
        for (obs_fn, obs_sam), (exp_fn, num_records, mapped_ids, unmapped_ids)\
                in zip(obs_sams, exp):
            self.assertEqual(str(obs_fn), exp_fn)
            with open(str(obs_sam)) as sam_f:
                obs_mapped_ids = [line.split('\t')[0]
                                  for line in sam_f
                                  if not line.startswith('@')]
            self.assertEqual(len(obs_mapped_ids), num_records)
            for e in mapped_ids:
                self.assertIn(e, obs_mapped_ids)
            for e in unmapped_ids:
                self.assertIn(e, obs_mapped_ids)

    def test_map_paired_alt_ceil_coefficient(self):
        obs_art, = self.plugin.methods['map_paired_reads'](
            self.demuxed_art, self.indexed_genome, ceil_coefficient=0.5)
        obs = obs_art.view(SAMFilesDirFmt)
        obs_sams = obs.sams.iter_views(SAMFormat)
        for _, obs_sam in obs_sams:
            with open(str(obs_sam)) as sam_f:
                self.assertIn('--n-ceil 0,0.5', sam_f.read())

    def test_sort_alignment_maps(self):
        obs_art, = self.plugin.methods['sort_alignment_maps'](
            self.paired_mapped_unsorted)
        obs = obs_art.view(SAMFilesDirFmt)
        exp_mapped_positions = \
            [1, 1, 192, 211, 402, 421, 612, 631, 823, 841]
        obs_sams = obs.sams.iter_views(SAMFormat)
        for _, obs_sam in obs_sams:
            with open(str(obs_sam)) as sam_f:
                obs_mapped_positions = \
                    [int(line.split('\t')[3])
                     for line in sam_f
                     if not line.startswith('@')]
                self.assertEqual(obs_mapped_positions,
                                 exp_mapped_positions)

    def test_remove_duplicates(self):
        obs_art, = self.plugin.methods['remove_duplicates'](
            self.sorted_alignment_maps)
        obs = obs_art.view(SAMFilesDirFmt)
        obs_sams = obs.sams.iter_views(SAMFormat)
        for _, obs_sam in obs_sams:
            with open(str(obs_sam)) as sam_f:
                obs_mapped_ids = [line.split('\t')[0]
                                  for line in sam_f
                                  if not line.startswith('@')]
            # one occurance of duplicate alignment is retained...
            self.assertIn('NB501727:157:HFWHJBGXF:2:22105:18312:6802',
                          obs_mapped_ids)
            # and the other isn't
            self.assertNotIn('NB501727:157:HFWHJBGXF:3:23610:2922:9731',
                             obs_mapped_ids)


if __name__ == '__main__':
    unittest.main()
