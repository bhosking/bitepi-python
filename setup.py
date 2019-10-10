from setuptools import setup, Extension


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

# Set __version__
exec(open('bitepi/version.py').read())
VERSION = __version__

setup(
    name='bitepi',
    version=VERSION,
    packages=[
        'bitepi',
    ],
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
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'License :: Other/Proprietary License',
        'Operating System :: POSIX',
        'Programming Language :: Python :: 3 :: Only',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
    python_requires='>=3'
)
