from distutils.core import setup
import os.path
import subprocess


VERSION = '0.1.1'


# Compile BitEpi
subprocess.run(
    args=[
        'g++',
        '-o',
        os.path.join('bitepi', 'BitEpi.o'),
        '-O3',
        os.path.join('source', 'BitEpi.cpp'),
        os.path.join('source', 'csvparser.c'),
        '-pthread',
    ],
    check=True,
)

# mark BitEpi.o executable
subprocess.run(
    args=[
        'chmod',
        '+x',
        os.path.join('bitepi', 'BitEpi.o'),
    ],
    check=True,
)

with open('README.md', 'r') as readme_file:
    long_description = readme_file.read()

setup(
    name='bitepi',
    version=VERSION,
    description=('Python wrapper for the BitEpi program, for finding epistasis'
                 'interactions'),
    long_description=long_description,
    author='Brendan Hosking',
    author_email='brendan.hosking@csiro.au',
    url='https://github.com/bhosking/bitepi-python',
    download_url=('https://github.com/bhosking/bitepi-python/archive/'
                  f'v{VERSION}.tar.gz'),
    keywords=['BIONINFORMATICS', 'EPISTASIS'],
    install_requires=[
        'numpy',
        'pandas',
    ],
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'License :: Other/Proprietary License',
        'Operating System :: POSIX',
        'Programming Language :: Python :: 3 :: Only',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
)

