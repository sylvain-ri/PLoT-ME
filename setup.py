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

__version__  = "0.8.5"
REQUIREMENTS = "requirements.txt"


with open("README.md", "r") as fh:
    long_description = fh.read()


def parse_requirements(path):
    list_pkg = []
    with open(path) as f:
        for line in f.readlines():
            if len(line) < 5 or "#" in line[:4]:
                continue
            else:
                list_pkg.append(line.replace(" ", "").replace("\n", ""))
    return list_pkg


setuptools.setup(
    name="PLoT-ME",
    version=__version__,
    author="Sylvain Riondet",
    author_email="sylvainrionder@gmail.com",
    description="Memory Reduction for Taxonomic Classifiers (pre-classifying long reads to clusters)",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/sylvain-ri/PLoT-ME",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3.7",
        "License :: OSI Approved :: MIT License",
        "Operating System :: POSIX :: Linux",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    python_requires='>=3.7',
    install_requires=parse_requirements(REQUIREMENTS),
    entry_points={
        'console_scripts': [
            'plot-me-parse = plot_me.parse_DB:arg_parser',
            'plot-me-classify = plot_me.classify:arg_parser',
        ],
    },
)
