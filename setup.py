from setuptools import setup, Extension, find_packages
import numpy as np

ext_modules = [
    Extension(
        "cpyrcolate.percolate_cy",
        sources=[
            "src/cpyrcolate/percolate_cy.c",
            "src/cpyrcolate/percolate_core.c",
        ],
        include_dirs=[np.get_include(), "src/cpyrcolate"],
        extra_compile_args=["-O3"],
    )
]

setup(
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    ext_modules=ext_modules,
    include_package_data=True,
)
