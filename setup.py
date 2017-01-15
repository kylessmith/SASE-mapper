from setuptools import setup, Extension
from Cython.Build import cythonize
import numpy as np

long_description = '''Software to identify regions of interest with a higher
                    than expected number of mutations in multiple samples'''

#extensions = [Extension('SASE_hunter/tpm',
#                        sources=['SASE_hunter/tpm/comb_pvals.pyx'],
#                        include_dirs=[np.get_include()]
#                        )]

setup(
    ext_modules = cythonize("SASE_mapper/tpm/tpm.pyx"),
    name="SASE_mapper",
    version="0.2.0",
    packages=["SASE_mapper"],
    author="Kyle S. Smith",
    license="MIT Licenses",
    description='Signatures of Accelerated Somatic Evolution mapper',
    install_requires=['numpy', 'scipy', 'pybedtools', 'bokeh', 'cython', 'matplotlib'],
    long_description=long_description,
    url="https://github.com/kylessmith/SASE-mapper",
    author_email="kyle.s.smith@ucdenver.edu"
)
