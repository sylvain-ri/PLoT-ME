#!/usr/bin/env python3
"""
Setup for pypi install
#############################################################################
Sylvain @ GIS / Biopolis / Singapore
Sylvain Jun-Zhe RIONDET <Riondet_Sylvain_from.tp@gis.a-star.edu.sg>
Started on 2019-12-11
PLoT-ME / Reads Binning Project
#############################################################################
"""

import setuptools

import plot_me

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="PLoT-ME",
    version=plot_me.__version__,
    author="Sylvain Riondet",
    author_email="sylvainrionder@gmail.com",
    description="PLoT-ME: Pre-Classification of Long reads for Memory Efficient Taxonomic Assignment",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/sylvain-ri/PLoT-ME",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3.7",
        "License :: MIT License",
        "Operating System :: Unix",
    ],
    python_requires='>=3.7',
    # todo: parse the requirements.txt
    install_requires=[],

    entry_points={
        'console_scripts': [
            'proclame-sm = sm_lib.core:proclamer',
        ],
    },
)
