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

from plot_me import __version__

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
    description="Pre-classification of Long-reads for Memory Efficient Taxonomic assignment",
    long_description=long_description,
    long_description_content_type="text/markdown",
    include_package_data=True,
    url="https://github.com/sylvain-ri/PLoT-ME",
    packages=setuptools.find_packages(exclude=("tests", "work_in_progress")),
    classifiers=[
        "Intended Audience :: Developers",
        "Intended Audience :: Healthcare Industry",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Operating System :: Unix",
        "Operating System :: POSIX :: Linux",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
        "Topic :: Database",
        "Topic :: Scientific/Engineering",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    python_requires='>=3.7',
    install_requires=parse_requirements(REQUIREMENTS),
    entry_points={
        'console_scripts': [
            'plot-me.preprocess = plot_me.parse_DB:arg_parser',
            'plot-me.classify = plot_me.classify:arg_parser',
        ],
    },
)
