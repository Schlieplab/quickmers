from setuptools import setup, Extension
import os
import numpy

## do this before running command: export CPATH=/opt/software/spack/spack/opt/spack/linux-rocky9-cascadelake/gcc-11.4.1/python-3.12.4-jxmctj6wbmmvwh77onnvucnn4wuxfr22/include/python3.12/:$CPATH
## and then pip install -e .

module = Extension(
    "quickmers._cbindings",
    sources=[
        "quickmers/_cbindings.c",
        "src/hamming.c",
        # add more C files here as you grow
    ],
    include_dirs=[
        numpy.get_include(),
        "quickmers/include"],
    extra_compile_args=["-O3"],
)

setup(
    name="quickmers",
    version="0.1.0",
    description="Fast k-mer utilities in C with Python bindings",
    author="Bioinformatics lab? change this", ##TODO: write correct author before publishing
    packages=["quickmers"],
    ext_modules=[module],
)