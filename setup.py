# https://packaging.python.org/tutorials/distributing-packages

from setuptools import setup
from codecs import open
import os

with open(os.path.join(os.path.abspath(os.path.dirname(__file__)), 'README.rst'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name='asap',
    use_scm_version={'write_to': 'asap/version.py'},
    setup_requires=['setuptools_scm'],
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
        'Programming Language :: Python :: 3 :: Only'
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.3',
        'Programming Language :: Python :: 3.4',
    ],
    entry_points = {
        'console_scripts': [
            'prepareJSONInput = asap.prepareJSONInput:main',
            'analyzeAmplicons = asap.analyzeAmplicons:main',
            'bamProcessor = asap.bamProcessor:main',
            'outputCombiner = asap.outputCombiner:main',
            'formatOutput = asap.formatOutput:main',
            'reformatXML = asap.reformatXML:main'
        ]
    },
    install_requires=[
        'numpy',
        'pysam<=0.13',
        'scikit-bio',
        'openpyxl',
        'lxml',
    ]
)
