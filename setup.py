# ----------------------------------------------------------------------------
# Copyright (c) 2020, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from setuptools import find_packages, setup

import versioneer


setup(
    name='q2-phylogenomics',
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    license='BSD-3-Clause',
    packages=find_packages(),
    author="Nicholas Bokulich",
    author_email="nbokulich@gmail.com",
    description="Phylogenomics toolkit and pipelines.",
    url="https://qiime2.org/",
    entry_points={
        'qiime2.plugins':
        ['q2-phylogenomics=q2_phylogenomics.plugin_setup:plugin']
    },
    package_data={
        'q2_phylogenomics.tests': ['data/*'],
        'q2_phylogenomics': ['citations.bib']
    },
    zip_safe=False,
)
