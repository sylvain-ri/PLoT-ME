# Cython accelerator
Most of the computing is done within a few methods of this project, so accelerating these reduces the overall runtime by a factor of 3 as of now.

https://cython.readthedocs.io/en/latest/

You will need to have Cython installed (tested with version 0.29).

If this doesn't work, PLoT-Me should automatically fall back to pure Python mode. There's a flag for both methods 
(preprocess and classify) to disable cython.

### Next steps to cythonize:
- binning of reads (import K-Means model)
- reading/writing of files

Building wheels is a bit tricky, but well explained there:
https://levelup.gitconnected.com/how-to-deploy-a-cython-package-to-pypi-8217a6581f09

In short: <br/>
`python3 -m pip install build auditwheel` <br/>
`python3 -m build` <br/>
`auditwheel repair <path>/PLoT_ME-<version>>-linux_x86_64.whl --plat manylinux2014_x86_64`<br/>


