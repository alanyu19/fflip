# -*- coding: utf-8 -*-

"""The setup script."""

from setuptools import setup, find_packages
import versioneer
import sys

# Custom commands for versioning
cmdclass = versioneer.get_cmdclass()

# Custom commands for building the docs
#try:
#    from sphinx.setup_command import BuildDoc
#    cmdclass['build_sphinx'] = BuildDoc
#except ImportError:
#    print("Could not import sphinx. build_sphinx command not installed.")


with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()


requirements = [
    'numpy>=1.6',
    'docutils>=0.18',
    'Click>=6.0,<8.0',
    'pandas',
    'nbsphinx', # sphinx documentation for jupyter notebooks
    'sphinx-click', # sphinx documentation for click
    'alabaster', # alabaster theme for documentation
    'packaging',
    'scipy>=1.0.0',
    # 'dask',
    'matplotlib',
    'plumbum',
    # 'dask-jobqueue>=0.7.0',
    # 'openmm>=7.3.0',
    # 'openmmtools',
    # 'mdtraj',
    # put package requirements here
]

if sys.version_info[0] < 3:
    requirements.append('matplotlib<3')
else:
    requirements.append('matplotlib')

setup_requirements = [
    'pytest-runner',
    'numpy>=1.20',
    # 'cython'
    #  put setup requirements (distutils extensions, etc.) here
]

test_requirements = [
    'pytest>=3.4',
    # put package test requirements here
]

setup(
    version=versioneer.get_version(),
    description="Force Field of Lipids Optimization Package",
    long_description=readme + '\n\n' + history,
    entry_points={
        'console_scripts': [
            'fflip=fflip.cli:entrypoint',
        ],
    },
    author="Yalun Yu",
    author_email='yalun.research@gmail.com',
    url='https://github.com/alanyu19/fflip',
    packages=find_packages(include=['fflip']),
    include_package_data=True,
    install_requires=requirements,
    license="GNU General Public License v3",
    zip_safe=False,
    keywords='fflip',
    name='fflip',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
    ],
    test_suite='tests',
    tests_require=test_requirements,
    setup_requires=setup_requirements,
    cmdclass=cmdclass,
)
