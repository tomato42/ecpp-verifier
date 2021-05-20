#!/usr/bin/env python

import io
import os

from setuptools import setup, find_packages

# Use README.md to set markdown long_description
directory = os.path.abspath(os.path.dirname(__file__))
readme_path = os.path.join(directory, "README.md")
with io.open(readme_path, encoding="utf-8") as read_file:
    long_description = read_file.read()

setup(
    name="ecpp",
    version="0.0.1",
    description="Tool to verify ECPP certificates",
    long_description=long_description,
    long_description_content_type="text/markdown",
    author="Hubert Kario",
    author_email="hkario@redhat.com",
    url="https://github.com/tomato42/ecpp-verifier",
    packages=find_packages("src"),
    scripts=["ecpp"],
    package_dir={"": "src"},
    include_package_data=True,
    exclude_package_data={"": [".gitignore"]},
    license="GPLv2.1",
    python_requires=">=2.7, !=3.0.*, !=3.1.*, !=3.2.*",
    classifiers=[
        "Programming Language :: Python",
        "Programming Language :: Python :: 2",
        "Programming Language :: Python :: 2.7",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.3",
        "Programming Language :: Python :: 3.4",
        "Programming Language :: Python :: 3.5",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
    ],
    install_requires=["six>=1.9.0", "ecdsa>=0.15"],
)
