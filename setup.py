from setuptools import setup, Extension, find_packages
from Cython.Build import cythonize
import numpy as np
import os

# Define the path to the C/Cython files
source_dir = os.path.join("src", "cpyrcolate")

ext_modules = [
    Extension(
        "cpyrcolate.percolate_cy",  # The importable compiled module
        sources=[
            os.path.join(source_dir, "percolate_cy.pyx"),
            os.path.join(source_dir, "percolate_core.c"),
        ],
        include_dirs=[np.get_include(), source_dir],
        extra_compile_args=["-O3"],  # Max C optimization
    )
]

setup(
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    ext_modules=cythonize(ext_modules, compiler_directives={"language_level": "3"}),
)
