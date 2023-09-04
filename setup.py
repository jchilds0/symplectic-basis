import os
from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
from setuptools.glob import glob

base_code = glob(os.path.join('snappy_src', 'kernel_code','*.c'))
unix_code = glob(os.path.join('snappy_src', 'unix_kit','*.c'))
addl_code = glob(os.path.join('snappy_src', 'addl_code', '*.c'))
symp_source_files = base_code + unix_code + addl_code

symp_source_files.append("cython/symp_basis.pyx")

symp_ext = Extension(
    name="symplectic_basis",
    sources=symp_source_files,
    include_dirs=["symp_src",
                  "snappy_src/addl_code",
                  "snappy_src/headers",
                  "snappy_src/kernel_code",
                  "snappy_src/real_type",
                  "snappy_src/unix_kit"],
    language="c"
)

setup(
    name="symplectic_basis",
    ext_modules=cythonize([symp_ext])
)