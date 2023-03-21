import os
import glob
import setuptools
from distutils.core import setup

with open("README.md", 'r') as readme:
    long_description = readme.read()

setup(
    name='vivarium-comets',
    version='0.0.1',
    packages=[
        'vivarium_comets',
        'vivarium_comets.processes',
        'vivarium_comets.composites',
        'vivarium_comets.experiments',
    ],
    author='Helen Scott',
    author_email='hscott@bu.edu',
    url='https://github.com/segrelab/vivarium-comets',
    license='MIT',
    entry_points={
        'console_scripts': []},
    short_description='a vivarium wrapper for COMETSpy',
    long_description=long_description,
    long_description_content_type='text/markdown',
    package_data={},
    include_package_data=True,
    install_requires=[
        'vivarium-core>=1.0.0',
        'pytest',
        # TODO: Add other dependencies.
    ],
)
