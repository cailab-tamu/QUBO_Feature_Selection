import pathlib
from setuptools import setup, find_packages

HERE = pathlib.Path(__file__).parent
README = (HERE / "README.md").read_text()
DESCRIPTION = "Quantum Feature Selection"

# Discover packages from the REPO ROOT. This finds the top-level `qfeatures`
# package (the folder with __init__.py). The old call find_packages(where=
# "qfeatures") searched *inside* qfeatures/ for sub-packages and returned [],
# so nothing was installed.
PACKAGES = find_packages(include=["qfeatures", "qfeatures.*"])
print(f"Packages found: {PACKAGES}")

# Match what the modules actually import. Keep versions unpinned so pip does
# not fight the conda env (environment.yml already pins the heavy ones).
INSTALL_REQUIRES = [
    "numpy",
    "pandas",
    "scipy",
    "scanpy",
    "joblib",
    "dwave-ocean-sdk",   # provides dwave.samplers, dwave.system, dimod
]
print(f"Install requires: {INSTALL_REQUIRES}")

KEYWORDS = [
    "quantum computing",
    "feature selection",
    "quantum annealing",
    "simulated annealing",
    "single-cell",
]

setup(
    name="qfeatures",
    version="0.1.0",
    long_description=README,
    long_description_content_type="text/markdown",
    description=DESCRIPTION,
    author="Selim Romero",
    url="https://github.com/cailab-tamu/QUBO_Feature_Selection",
    author_email="ssromerogon@tamu.edu",
    license="MIT",
    packages=PACKAGES,
    keywords=KEYWORDS,
    install_requires=INSTALL_REQUIRES,
    python_requires=">=3.9,<3.12",
)
