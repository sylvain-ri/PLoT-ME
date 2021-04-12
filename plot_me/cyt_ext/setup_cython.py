#! /usr/bin/env python3
# coding: utf-8
# defining NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

"""
Setup to build python files, with Cython

CMD: python3 v_array_setup.py build_ext --inplace

Not sure if this version is better than the recommended one from PyCharm:
https://www.jetbrains.com/help/pycharm/cython.html#get-started-cython
"""

from distutils.core import setup
from Cython.Build import cythonize

setup(
    name="PLoT-ME-cython",
    ext_modules=cythonize(['cyt_ext.pyx'],
                            language_level="3", ))
