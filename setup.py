import os
from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
from setuptools.glob import glob

base_code = glob(os.path.join('snappy_src/kernel_code', '*.c'))
symp_code = glob(os.path.join('symp_src', '*.c'))
symp_source_files = base_code + symp_code

symp_source_files.append("cython/symp_basis.pyx")

symp_ext = Extension(
    name="symplectic_basis",
    sources=symp_source_files,
    include_dirs=["symp_src",
                  "snappy_src/kernel_code",
                  "snappy_src/headers"],
    language="c"
)

setup(
    name="symplectic-basis",
    ext_modules=cythonize([symp_ext])
)