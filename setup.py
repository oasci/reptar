#!/usr/bin/env python

"""The setup script."""

from setuptools import setup, find_packages
import versioneer

with open('README.md') as readme_file:
    readme = readme_file.read()

requirements = [
    'numpy', 'exdir>=0.4.2', 'cclib>=1.7.0', 'scipy'
]

setup_requirements = ['versioneer']

test_requirements = requirements.append(['pytest', 'pytest-order'])

setup(
    author="Alex M. Maldonado",
    author_email='aalexmmaldonado@gmail.com',
    python_requires='>=3.7',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3'
    ],
    description="Implements reproducible data procedures",
    install_requires=requirements,
    extras_require={},
    license="MIT license",
    long_description=readme,
    include_package_data=True,
    package_data={'': ['definitions/*.yaml']},
    keywords='reptar',
    name='reptar',
    packages=find_packages(include=['reptar', 'reptar.*']),
    setup_requires=setup_requirements,
    test_suite='tests',
    tests_require=test_requirements,
    zip_safe=False,
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
)
