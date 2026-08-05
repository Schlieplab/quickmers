import sys
from setuptools import setup, Extension
import numpy

extra_args = ["-O3"]

if sys.platform == "win32":
    extra_args = ["/O2"]

ext = Extension(
    "quickmers._cbindings",
    sources=[
        "quickmers/_cbindings.c",
        "src/hamming.c",
        "src/levenshtein.c",
    ],
    include_dirs=[
        numpy.get_include(),
        "quickmers/include",
    ],
    extra_compile_args=extra_args,
)

setup(ext_modules=[ext])