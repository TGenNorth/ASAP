#!/usr/bin/env python3
# -*- coding: utf-8 -*-

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

readme = open('README.rst').read()
history = open('HISTORY.rst').read().replace('.. :changelog:', '')

requirements = [
    # TODO: put package requirements here
]

test_requirements = [
    # TODO: put package test requirements here
]

setup(
    name='ASAP',
    version='0.1',
    description='Amplicon Sequencing Analysis Pipeline',
    long_description=readme + '\n\n' + history,
    author='Darrin Lemmer',
    author_email='dlemmer@tgen.org',
    url='https://github.com/TGenNorth/ASAP',
    packages=[
        'asap'
    ],
    package_dir={
        'asap': 'asap',
    },
    include_package_data=True,
    install_requires=requirements,
    license="Academic and Research License",
    zip_safe=False,
    keywords='asap',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Developers',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.3',
        'Programming Language :: Python :: 3.4',
    ],
    #test_suite='tests',
    #tests_require=test_requirements,
    entry_points = {
        'console_scripts': [
            'prepareJSONInput = asap.prepareJSONInput:main',
            'analyzeAmplicons = asap.analyzeAmplicons:main',
            'bamProcessor = asap.bamProcessor:main',
            'outputCombiner = asap.outputCombiner:main',
            'formatOutput = asap.formatOutput:main'
        ]
    },
    scripts = [
    ]
)
