from setuptools import setup, find_packages

setup(
    name='qfeatures',  # The name of your package
    version='0.1.0',   # Version number
    description='A package for QUBO feature selection and related functions',
    author='Your Name',
    author_email='your_email@example.com',
    packages=find_packages(),  # Automatically find all packages in the directory
    install_requires=[
        'numpy',    # Add dependencies your package needs
        'pandas',
        'dwave-system',  # For D-Wave integration, if applicable
    ],
    python_requires='>=3.7',  # Specify Python version compatibility
)

