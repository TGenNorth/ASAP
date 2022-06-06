# https://packaging.python.org/tutorials/distributing-packages

from setuptools import setup
import codecs
from codecs import open
import os

with open(os.path.join(os.path.abspath(os.path.dirname(__file__)), 'README.rst'), encoding='utf-8') as f:
    long_description = f.read()

def read(rel_path):
    here = os.path.abspath(os.path.dirname(__file__))
    with codecs.open(os.path.join(here, rel_path), 'r') as fp:
        return fp.read()

def get_version(rel_path):
    for line in read(rel_path).splitlines():
        if line.startswith('__version__'):
            delim = '"' if '"' in line else "'"
            return line.split(delim)[1]
    else:
        raise RuntimeError("Unable to find version string.")

setup(
    name='asap',
    version=get_version("asap/__init__.py"),
    # use_scm_version={'write_to': 'asap/version.py'},
    # setup_requires=['setuptools_scm'],
    description='Amplicon Sequencing Analysis Pipeline',
    long_description=long_description,
    author='Darrin Lemmer',
    author_email='dlemmer@tgen.org',
    url='https://github.com/TGenNorth/ASAP',
    packages=[
        'asap'
    ],
    # include any data files it finds inside your package directories that are specified by your MANIFEST.in file
    include_package_data=True,
    license="Academic and Research License",
    zip_safe=False,
    keywords='bioinformatics',
    # Prevent pip from installing on other Python versions.
    # Requires setuptool >= 24.2.0 and pip >= 9.0.0
    python_requires='>=3',
    # https://pypi.python.org/pypi?%3Aaction=list_classifiers
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'Natural Language :: English',
        'Programming Language :: Python :: 3 :: Only',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.3',
        'Programming Language :: Python :: 3.4',
    ],
    entry_points = {
        'console_scripts': [
            'asap = asap.cmdParser:main',
        ]
    },
    install_requires=[
        'lxml',
        'numpy',
        'openpyxl',
        'pysam>0.13',
        'scikit-bio',
        'xmltodict',
    ]
)
