from setuptools import setup, Extension


VERSION = '0.1.7'

# the c++ extension module
bitepi_ext = Extension(
    name='bitepimodule',
    sources=[
        'bitepi/source/wrapper.cpp',
        'bitepi/source/csvparser.c',
    ],
    extra_compile_args=[
        '-O3',
        '-pthread',
    ],
)

with open('README.md', 'r') as readme_file:
    long_description = readme_file.read()

setup(
    name='bitepi',
    version=VERSION,
    description=('Python wrapper for the BitEpi program, for finding epistasis'
                 'interactions'),
    long_description=long_description,
    long_description_content_type='text/markdown',
    author='Brendan Hosking',
    author_email='brendan.hosking@csiro.au',
    url='https://github.com/bhosking/bitepi-python',
    download_url=('https://github.com/bhosking/bitepi-python/archive/'
                  'v{version}.tar.gz'.format(version=VERSION)),
    keywords=['BIONINFORMATICS', 'EPISTASIS'],
    install_requires=[
        'numpy',
        'pandas',
    ],
    ext_modules=[
        bitepi_ext,
    ],
    package_data={
        'bitepi': [
            'binary_source/*'
        ]
    },
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'License :: Other/Proprietary License',
        'Operating System :: POSIX',
        'Programming Language :: Python :: 3 :: Only',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
)
