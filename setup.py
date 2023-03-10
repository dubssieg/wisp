#!/usr/bin/env python3
from sys import version_info, stderr
from setuptools import setup, find_packages

NAME = 'wisp'
CURRENT_PYTHON = version_info[:2]
REQUIRED_PYTHON = (3, 10)

if CURRENT_PYTHON < REQUIRED_PYTHON:
    stderr.write(
        f"{NAME} requires Python 3.10 or higher and your current version is {CURRENT_PYTHON}.")
    exit(1)

with open("requirements.txt", "r", encoding='utf-8') as fh:
    requirements = [line.strip() for line in fh]

setup(
    name=NAME,
    version='0.1.0',
    description='XGBoost bacteria families identification',
    url='https://github.com/Tharos-ux/wisp',
    author='Tharos',
    author_email='dubois.siegfried@gmail.com',
    packages=['workspace'],  # find_packages(),
    package_data={'': ['parameters_files/params.json']},
    include_package_data=True,
    zip_safe=False,
    license="LICENSE",
    long_description=open("README.md", encoding='utf-8').read(),
    long_description_content_type='text/markdown',
    install_requires=requirements,
    entry_points={'console_scripts': [f'{NAME}=workspace.main:main']}
)
