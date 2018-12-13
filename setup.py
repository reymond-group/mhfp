import sys
import os
from setuptools import setup
from setuptools import Command
from setuptools.command.test import test as TestCommand
from datetime import datetime


with open('README.md', 'r') as fh:
    LONG_DESCRIPTION = fh.read()

NAME = 'mhfp'
VERSION = '1.4'
AUTHOR = 'Daniel Probst'
DESCRIPTION = 'Molecular MHFP fingerprints for cheminformatics applications'
URL = 'https://github.com/reymond-group/mhfp'
REQUIRED_PYTHON_VERSION = (3, 0)
PACKAGES = ['mhfp']
INSTALL_DEPENDENCIES = []
SETUP_DEPENDENCIES = [
]
TEST_DEPENDENCIES = [
    'pytest'
]
EXTRA_DEPENDENCIES = {
    'dev': [
        'pytest'
    ]
}

if sys.version_info < REQUIRED_PYTHON_VERSION:
    sys.exit('Python >= 3.0 is required. Your version:\n'+sys.version)


class PyTest(TestCommand):
    """
    Use pytest to run tests
    """
    user_options = [('pytest-args=', 'a', 'Arguments to pass into py.test')]

    def initialize_options(self):
        TestCommand.initialize_options(self)
        self.pytest_args = []

    def finalize_options(self):
        TestCommand.finalize_options(self)
        self.test_args = []
        self.test_suite = True

    def run_tests(self):
        import pytest
        errno = pytest.main(self.pytest_args)
        sys.exit(errno)


setup(
    name=NAME,
    description=DESCRIPTION,
    long_description=LONG_DESCRIPTION,
    long_description_content_type='text/markdown',
    url=URL,
    version=VERSION,
    author=AUTHOR,
    packages=PACKAGES,
    include_package_data=True,
    install_requires=INSTALL_DEPENDENCIES,
    setup_requires=SETUP_DEPENDENCIES,
    tests_require=TEST_DEPENDENCIES,
    extras_require=EXTRA_DEPENDENCIES,
    cmdclass={
        'test': PyTest
    }
)
