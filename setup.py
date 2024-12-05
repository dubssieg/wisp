#!/usr/bin/env python3
from sys import version_info, stderr
from setuptools import setup

NAME = 'wisp'
CURRENT_PYTHON = version_info[:2]
REQUIRED_PYTHON = (3, 10)

if CURRENT_PYTHON < REQUIRED_PYTHON:
    stderr.write(
        f"{NAME} requires Python 3.10 or higher and your current version is {CURRENT_PYTHON}.")
    exit(1)

setup(
    name=NAME,
    version='0.1.1',
    description='XGBoost bacteria families identification',
    url='https://github.com/Tharos-ux/wisp',
    author='Tharos',
    author_email='dubois.siegfried@gmail.com',
    packages=['workspace'],
    package_data={'': ['parameters_files/params.json']},
    include_package_data=True,
    zip_safe=False,
    license="MIT",
    long_description=open("README.md", encoding='utf-8').read(),
    long_description_content_type='text/markdown',
    install_requires=['tharos-pytools', 'xgboost', 'matplotlib', 'numpy',
                      'scipy', 'scikit-learn', 'seaborn', 'biopython', 'rich', 'treelib', 'plotly'],
    entry_points={'console_scripts': [f'{NAME}=workspace.main:main']}
)
