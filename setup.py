#!/usr/bin/env python3
"""
Setup for pypi install
#############################################################################
Sylvain @ GIS / Biopolis / Singapore
Sylvain Jun-Zhe RIONDET <Riondet_Sylvain_from.tp@gis.a-star.edu.sg>
Started on 2019-12-11
PLoT-ME / Reads Binning Project
#############################################################################

Setup to build python files, with Cython. For compilation testing:
CMD: python3 setup.py build_ext --inplace
For each new Cython file, with extension .pyx, add the header lines (see cyt_ext.pyx as example)
Cython FAQ : https://github.com/cython/cython/wiki/FAQ
If Pure Python is needed, then create .pxd file for each.py file
https://cython.readthedocs.io/en/latest/src/tutorial/pure.html

For Packaging, good tutorial for Cython package build:
https://levelup.gitconnected.com/how-to-deploy-a-cython-package-to-pypi-8217a6581f09
`python3 -m build`
`auditwheel repair dist/PLoT_ME-0.9.0-cp37-cp37m-linux_x86_64.whl --plat manylinux2014_x86_64`
add file ~/.pypirc : "[testpypi] \n username = __token__ \n password <token-value>"
`python3 -m twine upload --repository testpypi wheelhouse/*.whl`

For build and testing: python3 -m pip install -e .
"""
# Apparently setuptools is the news standard:
# https://stackoverflow.com/questions/32528560/
from setuptools import find_packages, setup, Extension
# setuptools MUST be imported first. https://stackoverflow.com/questions/21594925/
from Cython.Build import cythonize
from Cython.Distutils import build_ext
import numpy as np

# from plot_me import

REQUIREMENTS = "requirements.txt"


with open("README.md", "r") as fh:
    long_description = fh.read()


def parse_requirements(path):
    """ Get the list of dependencies from the requirements file """
    list_pkg = []
    with open(path) as f:
        for line in f.readlines():
            if len(line) < 5 or "#" in line[:4]:
                continue
            else:
                list_pkg.append(line.replace(" ", "").replace("\n", ""))
    return list_pkg


cython_module = cythonize([Extension("plot_me.cython_module.cyt_ext", ['plot_me/cython_module/cyt_ext.pyx'],
                                    extra_compile_args=['-fopenmp'],
                                    extra_link_args=['-fopenmp'],
                                    include_dirs=["."],), ],
                          language_level="3")
# Can use Extension()

setup(
    name="PLoT-ME",
    version="0.9.0",
    author="Sylvain Riondet",
    author_email="sylvainrionder@gmail.com",
    description="Pre-classification of Long-reads for Memory Efficient Taxonomic assignment",
    long_description=long_description,
    long_description_content_type="text/markdown",
    include_package_data=True,
    url="https://github.com/sylvain-ri/PLoT-ME",
    packages=find_packages(exclude=("tests", "work_in_progress", "plot_me/cython_module/work_in_progress")),
    include_dirs=np.get_include(),
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
            # todo: Direct access to kmer counter
        ],
    },
    ext_modules=cython_module,
    cmdclass={'build_ext': build_ext},
)
