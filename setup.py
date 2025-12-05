from setuptools import setup, Extension
import numpy

## do this before running command: export CPATH=/opt/software/spack/spack/opt/spack/linux-rocky9-cascadelake/gcc-11.4.1/python-3.12.4-jxmctj6wbmmvwh77onnvucnn4wuxfr22/include/python3.12/:$CPATH
## and then pip install -e .

with open("README.md", "r") as f:
    long_description = f.read()

module = Extension(
    "quickmers._cbindings",
    sources=[
        "quickmers/_cbindings.c",
        "src/hamming.c",
        "src/levenshtein.c",
        # add more C files here as you grow
    ],
    include_dirs=[
        numpy.get_include(),
        "quickmers/include"],
    extra_compile_args=["-O3", "-mavx2", "-march=native"]
)

setup(
    name="quickmers",
    version="0.1.0",
    description="Fast k-mer utilities in C with Python bindings",
    author="Bioinformatics lab? change this", ##TODO: write correct author before publishing
    author_email="TODO", ##TODO: write correct email before publishing
    long_description=long_description,
    long_description_content_type="text/markdown",
    packages=["quickmers"],
    ext_modules=[module],
    url='https://gitea.schlieplab.org/Research/QuickMers' ##TODO: change this to github link
)