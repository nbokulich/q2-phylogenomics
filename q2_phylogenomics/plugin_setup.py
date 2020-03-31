# ----------------------------------------------------------------------------
# Copyright (c) 2020, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from qiime2.plugin import (
    Choices,
    Plugin,
    Citations,
    Range,
    Int,
    Str,
    List,
    Bool,
)
from q2_types.feature_data import FeatureData, Sequence
from q2_types.sample_data import SampleData
from q2_types.per_sample_sequences import (
    SequencesWithQuality,
    PairedEndSequencesWithQuality,
)

import q2_phylogenomics
import q2_phylogenomics._prinseq
import q2_phylogenomics._filter
from q2_types.bowtie2 import Bowtie2Index


citations = Citations.load('citations.bib', package='q2_phylogenomics')

plugin = Plugin(
    name='phylogenomics',
    version=q2_phylogenomics.__version__,
    website='https://github.com/qiime2/q2-phylogenomics',
    package='q2_phylogenomics',
    description='A QIIME 2 plugin for phylogenomics analyses.',
    short_description='A QIIME 2 plugin for phylogenomics analyses.',
)

prinseq_input = {'demultiplexed_sequences': 'The sequences to be trimmed.'}
prinseq_output = {'trimmed_sequences': 'The resulting trimmed sequences.'}

prinseq_parameters = {
    'trim_qual_right': Int % Range(1, None),
    'trim_qual_type': Str % Choices(['min', 'mean', 'max', 'sum']),
    'trim_qual_window': Int % Range(1, None),
    'min_qual_mean': Int % Range(1, None),
    'min_len': Int % Range(1, None),
    'lc_method': Str % Choices(['dust', 'entropy']),
    'lc_threshold': Int % Range(0, 100),
    'derep': List[Str % Choices(list('12345'))]}

prinseq_parameter_descriptions = {
    'trim_qual_right': 'Trim sequence by quality score from the 3\'-end with '
                       'this threshold score.',
    'trim_qual_type': 'Type of quality score calculation to use. Allowed '
                      'options are min, mean, max and sum.',
    'trim_qual_window': 'The sliding window size used to calculate quality '
                        'score by type. To stop at the first base that fails '
                        'the rule defined, use a window size of 1.',
    'min_qual_mean': 'Filter sequence with quality score mean below '
                     'min_qual_mean.',
    'min_len': 'Filter sequence shorter than min_len.',
    'lc_method': 'Method to filter low complexity sequences.',
    'lc_threshold': 'The threshold value used to filter sequences by sequence '
                    'complexity. The dust method uses this as maximum allowed '
                    'score and the entropy method as minimum allowed value.',
    'derep': 'Type of duplicates to filter. Use integers for multiple '
             'selections (e.g. 124 to use type 1, 2 and 4). The order does '
             'not matter. Option 2 and 3 will set 1 and option 5 will set 4 '
             'as these are subsets of the other option.\n\n1 (exact '
             'duplicate), 2 (5\' duplicate), 3 (3\' duplicate), 4 (reverse '
             'complement exact duplicate), 5 (reverse complement 5\'/3\' '
             'duplicate).'
}

plugin.methods.register_function(
    function=q2_phylogenomics._prinseq.prinseq_single,
    inputs={'demultiplexed_sequences': SampleData[SequencesWithQuality]},
    parameters=prinseq_parameters,
    outputs=[('trimmed_sequences', SampleData[SequencesWithQuality])],
    input_descriptions=prinseq_input,
    parameter_descriptions=prinseq_parameter_descriptions,
    output_descriptions=prinseq_output,
    name='Filter and trim demultiplexed single-end sequences with PRINSEQ.',
    description='Filter and trim demultiplexed single-end FASTQ sequences '
                'based on quality scores using PRINSEQ-lite.',
    citations=[citations['schmieder_prinseq']]
)

plugin.methods.register_function(
    function=q2_phylogenomics._prinseq.prinseq_paired,
    inputs={
        'demultiplexed_sequences': SampleData[PairedEndSequencesWithQuality]},
    parameters=prinseq_parameters,
    outputs=[('trimmed_sequences', SampleData[PairedEndSequencesWithQuality])],
    input_descriptions=prinseq_input,
    parameter_descriptions=prinseq_parameter_descriptions,
    output_descriptions=prinseq_output,
    name='Filter and trim demultiplexed paired-end sequences with PRINSEQ.',
    description='Filter and trim demultiplexed paired-end FASTQ sequences '
                'based on quality scores using PRINSEQ-lite.',
    citations=[citations['schmieder_prinseq']]
)

filter_input = {'demultiplexed_sequences': 'The sequences to be trimmed.',
                'database': 'Bowtie2 indexed database.'}
filter_output = {'filtered_sequences': 'The resulting filtered sequences.'}

filter_parameters = {
    'n_threads': Int % Range(1, None),
    'mode': Str % Choices(['local', 'global']),
    'sensitivity': Str % Choices([
        'very-fast', 'fast', 'sensitive', 'very-sensitive']),
    'ref_gap_open_penalty': Int % Range(1, None),
    'ref_gap_ext_penalty': Int % Range(1, None),
    'exclude_seqs': Bool,
}

filter_parameter_descriptions = {
    'n_threads': 'Number of alignment threads to launch.',
    'mode': 'Bowtie2 alignment settings. See bowtie2 manual for more details.',
    'sensitivity': 'Bowtie2 alignment sensitivity. See bowtie2 manual for '
                   'details.',
    'ref_gap_open_penalty': 'Reference gap open penalty.',
    'ref_gap_ext_penalty': 'Reference gap extend penalty.',
    'exclude_seqs': 'Exclude sequences that align to reference. Set this '
                    'option to False to exclude sequences that do not align '
                    'to the reference database.'
}

filter_citations = [citations['langmead2012fast'],
                    citations['heng2009samtools']]

filter_description = (
    'Filter out (or keep) sequences that align to reference database, using '
    'bowtie2 and samtools. This method can be used to filter out human DNA '
    'sequences and other contaminant in any FASTQ sequence data (e.g., '
    'shotgun genome or amplicon sequence data), or alternatively (when '
    'exclude_seqs is False) to only keep sequences that do align to the '
    'reference.')

plugin.methods.register_function(
    function=q2_phylogenomics._filter.filter_single,
    inputs={'demultiplexed_sequences': SampleData[SequencesWithQuality],
            'database': Bowtie2Index},
    parameters=filter_parameters,
    outputs=[('filtered_sequences', SampleData[SequencesWithQuality])],
    input_descriptions=filter_input,
    parameter_descriptions=filter_parameter_descriptions,
    output_descriptions=filter_output,
    name='Filter single-end sequences by alignment to reference database.',
    description=filter_description,
    citations=filter_citations
)

plugin.methods.register_function(
    function=q2_phylogenomics._filter.filter_paired,
    inputs={
        'demultiplexed_sequences': SampleData[PairedEndSequencesWithQuality],
        'database': Bowtie2Index},
    parameters=filter_parameters,
    outputs=[
        ('filtered_sequences', SampleData[PairedEndSequencesWithQuality])],
    input_descriptions=filter_input,
    parameter_descriptions=filter_parameter_descriptions,
    output_descriptions=filter_output,
    name='Filter paired-end sequences by alignment to reference database.',
    description=filter_description,
    citations=filter_citations
)

plugin.methods.register_function(
    function=q2_phylogenomics._filter.bowtie2_build,
    inputs={'sequences': FeatureData[Sequence]},
    parameters={'n_threads': Int % Range(1, None)},
    outputs=[('database', Bowtie2Index)],
    input_descriptions={
        'sequences': 'Reference sequences used to build bowtie2 index.'},
    parameter_descriptions={'n_threads': 'Number of threads to launch'},
    output_descriptions={'database': 'Bowtie2 index.'},
    name='Build bowtie2 index from reference sequences.',
    description='Build bowtie2 index from reference sequences.',
    citations=[citations['langmead2012fast']]
)
