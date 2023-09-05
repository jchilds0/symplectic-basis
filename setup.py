import os
from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
from setuptools.glob import glob

base_code = glob(os.path.join('snappy_src/kernel_code', '*.c'))
symp_code = glob(os.path.join('symp_src', '*.c'))
symp_source_files = base_code + symp_code

symp_source_files.append("cython/symp_basis.pyx")

base_head = glob(os.path.join('snappy_src/headers', '*.h'))
symp_head = glob(os.path.join('symp_src', '*.h'))

symp_header_files = base_head + symp_head

symp_ext = Extension(
    name="symplectic_basis",
    sources=symp_source_files,
    library_dirs=["symp_src",
                  "snappy_src/kernel_code",
                  "snappy_src/headers"],
    include_dirs=["symp_src",
                  "snappy_src/kernel_code",
                  "snappy_src/headers"],
    language="c"
)

setup(
    name="symplectic-basis",
    headers=symp_header_files,
    ext_modules=cythonize(symp_ext, compiler_directives={'language_level': "3"})
)