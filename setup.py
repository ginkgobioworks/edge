#!/usr/bin/env python

from setuptools import setup, find_packages

version = "0.2"

with open("requirements/core.txt", "r") as f:
    install_requires = [x.strip() for x in f.readlines() if not x.strip().startswith("#")]

setup(name="edge",
      version=version,
      author="Ginkgo Bioworks",
      author_email="team@ginkgobioworks.com",
      description="Genome Engineering Tool",
      license="MIT",
      namespace_packages=["src"],
      packages=find_packages(),
      zip_safe=False,
      install_requires=install_requires)
