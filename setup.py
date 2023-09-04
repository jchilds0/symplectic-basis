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
    name="symplectic-basis",
    version='0.0.1',
    license='GPLv3+',
    description='A symplectic basis for triangulated 3-manifolds, accepts manifolds from SnapPy',
    author='Josh Childs',
    url='https://github.com/jchilds0/symplectic-basis',
    keywords=["SnapPy", "3-manifolds", "Hyperbolic Geometry"],
    download_url="https://github.com/jchilds0/symplectic-basis/archive/refs/tags/v0.0.1.tar.gz",
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Mathematics',
        'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',
        'Programming Language :: Python :: 3',
    ],
    ext_modules=cythonize([symp_ext])
)