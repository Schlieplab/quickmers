import platform
import subprocess
from setuptools import setup, Extension
import numpy

def cpu_supports_avx2():
    if platform.system() != "Linux" and platform.system() != "Darwin":
        return False
    try:
        output = subprocess.check_output("lscpu", shell=True, text=True)
        return "avx2" in output.lower()
    except Exception:
        return False

extra_args = ["-O3", "-march=native"]
if cpu_supports_avx2():
    extra_args.append("-mavx2")

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
    extra_compile_args=extra_args
)

setup(ext_modules=[ext])
