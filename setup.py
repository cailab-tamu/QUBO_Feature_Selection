import pathlib
from setuptools import setup, find_packages

# Ensure the path to the 'qfeatures' folder is correctly added
HERE = pathlib.Path(__file__).parent
README = (HERE / "README.md").read_text()
DESCRIPTION = "Quantum Feature Selection"

# Find packages inside the 'qfeatures' folder
PACKAGES = find_packages(where="qfeatures")

# Verbose output: print the list of packages found
print(f"Packages found: {PACKAGES}")

INSTALL_REQUIRES = [
    'numpy',    # Add dependencies your package needs
    'pandas',
    'dwave-system',  # For D-Wave integration, if applicable
]

# Print the install_requires for visibility
print(f"Install requires: {INSTALL_REQUIRES}")

KEYWORDS = [
    "quantum computing",
    "feature selection",
    "quantum annealing",
    "simulated annealing",
    "single-cell"
]

setup(
    name='qfeatures',  # Updated package name
    version='0.1.0',   # Version number
    long_description=README,
    description=DESCRIPTION,
    author='Selim Romero',
    url='https://github.com/cailab-tamu/QUBO_Feature_Selection',
    author_email='ssromerogon@tamu.edu',
    license="MIT",
    packages=PACKAGES,  # Look for packages inside the 'qfeatures' folder
    keywords=KEYWORDS,
    package_dir={'qfeatures': 'qfeatures'},  # Adjust to match your structure
    install_requires=INSTALL_REQUIRES,  # Add dependencies
    python_requires='>=3.9,<3.10',  # Specify Python version compatibility for any 3.9.x version
)
