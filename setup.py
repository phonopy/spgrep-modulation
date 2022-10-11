#!/usr/bin/env python
# https://github.com/kennethreitz/setup.py/blob/master/setup.py

import os

from setuptools import find_packages, setup

# Package meta-data.
NAME = "spgrep_modulation"
DESCRIPTION = "Collective atomic modulation analysis with irreducible space-group representation"
URL = "https://github.com/lan496/spgrep-modulation"
AUTHOR = "Kohei Shinohara"
EMAIL = "kshinohara0508@gmail.com"
REQUIRES_PYTHON = ">=3.8.0"

# What packages are required for this module to be executed?
REQUIRED = [
    "setuptools",
    "setuptools_scm",
    "wheel",
    "typing_extensions",
    "numpy>=1.20.1",
    "spglib>=1.16.5",
    "phonopy>=2.15.1",
    "seekpath",
    "spgrep>=0.2.11",  # https://github.com/spglib/spgrep
    "hsnf>=0.3.15",  # https://github.com/lan496/hsnf
]

VIS_REQUIRED = [
    "ipykernel",
    "notebook>=4.2",
    # nglview does not work with ipywidgets>=8, https://github.com/nglviewer/nglview/issues/1032
    "ipywidgets>=7.0, <8",
    "nglview>=3.0.3",
    "ase>=3.22.1",
    "pymatgen>=2022.8.23",
]

# What packages are optional?
EXTRAS = {
    "vis": VIS_REQUIRED,
    "dev": [
        "pytest",
        "pytest-cov",
        "pre-commit",
        "black",
        "mypy",
        "flake8",
        "pyupgrade",
        "ipython",
        "click",
        "networkx",
        *VIS_REQUIRED,
    ],
    "docs": [
        "sphinx",
        "sphinx-autobuild",
        "nbsphinx",
        "sphinxcontrib-bibtex",
        # "sphinx-proof",  # Seems to be incompatible with nbsphinx
        "myst-parser",
        "sphinx-book-theme",
    ],
}

# The rest you shouldn't have to touch too much :)
# ------------------------------------------------
# Except, perhaps the License and Trove Classifiers!
# If you do change the License, remember to change the Trove Classifier for that!

here = os.path.abspath(os.path.dirname(__file__))

# Import the README and use it as the long-description.
# Note: this will only work if 'README.md' is present in your MANIFEST.in file!
try:
    with open(os.path.join(here, "README.md"), encoding="utf-8") as f:
        long_description = "\n" + f.read()
except FileNotFoundError:
    long_description = DESCRIPTION


# Where the magic happens:
setup(
    name=NAME,
    description=DESCRIPTION,
    long_description=long_description,
    long_description_content_type="text/markdown",
    author=AUTHOR,
    author_email=EMAIL,
    python_requires=REQUIRES_PYTHON,
    url=URL,
    package_dir={"": "src"},
    packages=find_packages(where="src", include=["spgrep_modulation"]),
    package_data={},
    # If your package is a single module, use this instead of 'packages':
    # py_modules=['mypackage'],
    # entry_points={
    #     'console_scripts': ['mycli=mymodule:cli'],
    # },
    # numpy: https://github.com/numpy/numpy/issues/2434
    setup_requires=["setuptools_scm", "numpy"],
    install_requires=REQUIRED,
    extras_require=EXTRAS,
    include_package_data=True,
    license="BSD",
    test_suite="tests",
    zip_safe=False,
    use_scm_version=True,
    classifiers=[
        # Trove classifiers
        # Full list: https://pypi.python.org/pypi?%3Aaction=list_classifiers
        "Development Status :: 2 - Pre-Alpha",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: BSD License",
        "Operating System :: OS Independent",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Topic :: Scientific/Engineering :: Physics",
    ],
)
