import io
import os
from setuptools import setup, find_packages

# Package meta-data.
NAME = 'bacdiving'
DESCRIPTION = 'Bacdiving accesses the Bacterial Diversity Metadatabase BacDive and provides various visualization options.'
URL = 'https://github.com/mBiocoder/Bacdiving'
EMAIL = 'M.Arunkumar@campus.lmu.de'
AUTHOR = 'Mahima Arunkumar'
REQUIRES_PYTHON = '>=3.8.0'
VERSION = "1.2.3"

# What packages are required for this module to be executed?
REQUIRED = [
    'alive_progress>=2.4.1',
    'anytree>=2.8.0',
    'bacdive>=0.2',
    'bokeh>=2.4.3',
    'ete3>=3.1.2',
    'matplotlib==3.6.0',
    'numpy>=1.23.0',
    'pandas>=1.5.0',
    'seaborn',
    'setuptools>=65.5.0',
    'scipy>=1.9.2',
    'toyplot>=1.0.3',
    'toytree>=2.0.1',
    'wheel>=0.34.1',
    'worldmap>=0.1.6']



long_description = 'README.md'



# Where the magic happens:
setup(
    name=NAME,
    version=VERSION,
    description=DESCRIPTION,
    long_description=long_description,
    long_description_content_type='text/markdown',
    author=AUTHOR,
    author_email=EMAIL,
    python_requires=REQUIRES_PYTHON,
    url=URL,
    install_requires=REQUIRED,
    include_package_data=True,
    license='MIT',
    keywords=["BacDive", "bacteria", "phenotype information"]
)
