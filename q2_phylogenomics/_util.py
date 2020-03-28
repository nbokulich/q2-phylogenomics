# ----------------------------------------------------------------------------
# Copyright (c) 2020, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------


import subprocess
import gzip
import shutil


def run_command(cmd, verbose=True):
    print('Running external command line application. This may print '
          'messages to stdout and/or stderr.')
    print('The commands to be run are below. These commands cannot '
          'be manually re-run as they will depend on temporary files that '
          'no longer exist.')
    print('\nCommand:', end=' ')
    print(' '.join(cmd), end='\n\n')
    subprocess.run(cmd, check=True)


def _gzip_compress(input_fp, output_fp):
    with open(input_fp, 'rb') as temp_in:
        with gzip.open(output_fp, 'wb') as temp_out:
            shutil.copyfileobj(temp_in, temp_out)


def _gzip_decompress(input_fp, output_fp):
    with gzip.open(input_fp, 'rt') as temp_in:
        with open(output_fp, 'w') as temp_out:
            shutil.copyfileobj(temp_in, temp_out)
