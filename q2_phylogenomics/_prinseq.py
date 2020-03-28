# ----------------------------------------------------------------------------
# Copyright (c) 2020, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import shutil
import tempfile
import pandas as pd

from q2_types.per_sample_sequences import (
    CasavaOneEightSingleLanePerSampleDirFmt,
    SingleLanePerSampleSingleEndFastqDirFmt,
    SingleLanePerSamplePairedEndFastqDirFmt,
)

from ._util import run_command, _gzip_compress, _gzip_decompress


_prinseq_defaults = {
    'trim_qual_right': 30,
    'trim_qual_type': 'min',
    'trim_qual_window': 5,
    'min_qual_mean': 20,
    'min_len': 70,
    'lc_method': 'dust',
    'lc_threshold': 3,
    'derep': ['1', '4']
}


def _run_prinseq(
        f_read, r_read, trimmed_seqs,
        trim_qual_right=_prinseq_defaults['trim_qual_right'],
        trim_qual_type=_prinseq_defaults['trim_qual_type'],
        trim_qual_window=_prinseq_defaults['trim_qual_window'],
        min_qual_mean=_prinseq_defaults['min_qual_mean'],
        min_len=_prinseq_defaults['min_len'],
        lc_method=_prinseq_defaults['lc_method'],
        lc_threshold=_prinseq_defaults['lc_threshold'],
        derep=_prinseq_defaults['derep'],
        ):
    derep = ''.join(derep)
    # prinseq-lite only accepts unzipped fastq
    temp_dir = tempfile.mkdtemp(prefix='a-place-to-put-unzipped-fastqs-')
    f_out = '{0}/{1}.fastq'.format(temp_dir, os.path.basename(f_read))
    _gzip_decompress(f_read, f_out)
    if r_read is not None:
        r_out = '{0}/{1}.fastq'.format(temp_dir, os.path.basename(r_read))
        _gzip_decompress(r_read, r_out)

    outname = temp_dir + '/outfile'

    cmd = [
        'prinseq-lite.pl',
        '-trim_qual_right', str(trim_qual_right),
        '-trim_qual_type', str(trim_qual_type),
        '-trim_qual_window', str(trim_qual_window),
        '-min_qual_mean', str(min_qual_mean),
        '-min_len', str(min_len),
        '-lc_method', str(lc_method),
        '-lc_threshold', str(lc_threshold),
        '-derep', str(derep),
        '-out_good', outname,
        '-out_bad', 'null',
        '-fastq', f_out,
    ]

    if r_read is not None:
        cmd += ['-fastq2', r_out]

    run_command(cmd)

    # copy prinseq output to its new home
    # prinseq has its own output path naming scheme, so rename to keep Q2 happy
    if r_read is not None:
        r_read = str(trimmed_seqs.path / os.path.basename(r_read))
        _gzip_compress(outname + '_2.fastq', r_read)
        # if using paired-end data, adjust the trimmed forward read filepath.
        # For details see prinseq-lite manual -out_good option.
        outname += '_1'
    f_read = str(trimmed_seqs.path / os.path.basename(f_read))
    _gzip_compress(outname + '.fastq', f_read)

    shutil.rmtree(temp_dir)


def prinseq_single(
        demultiplexed_sequences: SingleLanePerSampleSingleEndFastqDirFmt,
        trim_qual_right: int = _prinseq_defaults['trim_qual_right'],
        trim_qual_type: str = _prinseq_defaults['trim_qual_type'],
        trim_qual_window: int = _prinseq_defaults['trim_qual_window'],
        min_qual_mean: int = _prinseq_defaults['min_qual_mean'],
        min_len: int = _prinseq_defaults['min_len'],
        lc_method: str = _prinseq_defaults['lc_method'],
        lc_threshold: int = _prinseq_defaults['lc_threshold'],
        derep: str = _prinseq_defaults['derep']) -> \
            CasavaOneEightSingleLanePerSampleDirFmt:
    trimmed_sequences = CasavaOneEightSingleLanePerSampleDirFmt()
    df = demultiplexed_sequences.manifest.view(pd.DataFrame)
    for _, fwd in df.itertuples():
        _run_prinseq(fwd, None, trimmed_sequences, trim_qual_right,
                     trim_qual_type, trim_qual_window, min_qual_mean, min_len,
                     lc_method, lc_threshold, derep)
    return trimmed_sequences


def prinseq_paired(
        demultiplexed_sequences: SingleLanePerSamplePairedEndFastqDirFmt,
        trim_qual_right: int = _prinseq_defaults['trim_qual_right'],
        trim_qual_type: str = _prinseq_defaults['trim_qual_type'],
        trim_qual_window: int = _prinseq_defaults['trim_qual_window'],
        min_qual_mean: int = _prinseq_defaults['min_qual_mean'],
        min_len: int = _prinseq_defaults['min_len'],
        lc_method: str = _prinseq_defaults['lc_method'],
        lc_threshold: int = _prinseq_defaults['lc_threshold'],
        derep: str = _prinseq_defaults['derep']) -> \
            CasavaOneEightSingleLanePerSampleDirFmt:
    trimmed_sequences = CasavaOneEightSingleLanePerSampleDirFmt()
    df = demultiplexed_sequences.manifest.view(pd.DataFrame)
    for _, fwd, rev in df.itertuples():
        _run_prinseq(fwd, rev, trimmed_sequences, trim_qual_right,
                     trim_qual_type, trim_qual_window, min_qual_mean, min_len,
                     lc_method, lc_threshold, derep)
    return trimmed_sequences
