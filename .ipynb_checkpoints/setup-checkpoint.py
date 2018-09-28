#!/usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import setup, find_packages

with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

requirements = [
    "pyensembl"
]

test_requirements = [
    "py.test"
]

setup(
    name='gfeat',
    version='0.1.0',
    description="Python genomic features extractor from raw files.",
    long_description=readme + '\n\n' + history,
    author="Veronika Kotova",
    author_email='kotova@in.tum.de',
    url='https://github.com/gagneurlab/gfeat',
    packages=find_packages(),
    include_package_data=True,
    install_requires=requirements,
    license="MIT license",
    zip_safe=False,
    keywords='gfeat',
    extras_require={
        "develop": ["bumpversion",
                    "wheel",
                    "jedi",
                    "epc",
                    "pytest",
                    "pytest-pep8",
                    "pytest-cov"],
    },
    classifiers=[
        'Programming Language :: Python :: 3.5',
    ],
    test_suite='tests',
    tests_require=test_requirements
)
