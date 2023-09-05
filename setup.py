import os
from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
from setuptools.glob import glob

base_code = glob(os.path.join('snappea_src', 'kernel_code','*.c'))
unix_code = glob(os.path.join('snappea_src', 'unix_kit','*.c'))
addl_code = glob(os.path.join('snappea_src', 'addl_code', '*.c'))
symp_code = glob(os.path.join('symp_src', '*.c'))
symp_source_files = base_code + unix_code + addl_code + symp_code

symp_source_files.append("cython/symp_basis.pyx")

symp_ext = Extension(
    name="symplectic_basis",
    sources=symp_source_files,
    include_dirs=["symp_src",
                  "snappea_src/addl_code",
                  "snappea_src/headers",
                  "snappea_src/kernel_code",
                  "snappea_src/real_type",
                  "snappea_src/unix_kit"],
    language="c"
)

setup(
    name="symplectic-basis",
    ext_modules=cythonize([symp_ext])
)