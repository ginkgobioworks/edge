#!/usr/bin/env python

import glob
import os
import pip
from setuptools import setup
from setuptools.command.install import install as InstallCommand

PACKAGE_ROOT = os.path.relpath(os.path.abspath(os.path.split(__file__)[0]), os.getcwd())

class EdgeInstallCommand(InstallCommand):
    deployment_options = ("production", "development")
    deployment_options_txt = "deployment can be %r" % str.join(', ', deployment_options)
    user_options = [
        ("deployment=", "d", deployment_options_txt),
    ] + InstallCommand.user_options

    def initialize_options(self):
        InstallCommand.initialize_options(self)
        self.deployment = "production"

    def finalize_options(self):
        for deployment_option in self.deployment_options:
            if deployment_option.startswith(self.deployment.lower()):
                self.deployment = deployment_option
                break
        assert self.deployment in self.deployment_options, self.deployment_options_txt
        self.populate_requirements()
        InstallCommand.finalize_options(self)

    def populate_requirements(self):
        reqfn = os.path.join(PACKAGE_ROOT, "requirements", "%s.txt" % self.deployment)
        requirements = pip.req.parse_requirements(reqfn, session=pip.download.PipSession())
        self.distribution.install_requires = [str(req.req) for req in requirements]

def gather_scripts(config):
    scripts = glob.glob(os.path.join(PACKAGE_ROOT, "scripts", "edge_*"))
    config["scripts"] = scripts

def setup_config():
    config = {
        "name": "edge",
        "version": "0.3",
        "author": "Ginkgo Bioworks",
        "author_email": "team@ginkgobioworks.com",
        "description": "Genome Engineering Tool",
        "license": "MIT",
        "packages": ["edge", "edge.server"],
        "package_dir": {"edge": "src"},
        "include_package_data": True,
        "zip_safe": False,
        "cmdclass": {"install": EdgeInstallCommand},
    }
    gather_scripts(config)
    return config

if __name__ == "__main__":
    edge_setup_config = setup_config()
    setup(**edge_setup_config)
