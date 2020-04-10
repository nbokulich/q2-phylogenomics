# ----------------------------------------------------------------------------
# Copyright (c) 2020, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------


from qiime2.plugin import SemanticType

from q2_types.sample_data import SampleData


AlignmentMap = SemanticType("AlignmentMap",
                            variant_of=SampleData.field['type'])
PileUp = SemanticType("PileUp", variant_of=SampleData.field['type'])
ConsensusSequences = SemanticType("ConsensusSequences",
                                  variant_of=SampleData.field['type'])
ReferenceSequence = SemanticType("ReferenceSequence")
