#!/usr/bin/env python

import glob
import os
import pip
import setuptools
from setuptools.command.install import install as InstallCommand

SourceDir = os.path.split(os.path.abspath(__file__))[0]

class BowerBuildCommand(setuptools.Command):
    description = "run bower commands."
    user_options = [
        ('bower-command=', 'c',
         'Bower command to run. Defaults to \'install\'.'),
        ('force-latest', 'F', 'Force latest version on conflict.'),
        ('production', 'p', 'Do not install project devDependencies.'),
    ]

    boolean_options = ['production', 'force-latest']

    def initialize_options(self):
        self.force_latest = False
        self.production = False
        self.bower_command = 'install'

    def finalize_options(self):
        pass

    def run(self):
        cmd = ['bower', self.bower_command]
        if self.force_latest:
            cmd.append('-F')
        if self.production:
            cmd.append('-p')
        self.spawn(cmd)


class EdgeInstallCommand(InstallCommand):
    deployment_options = ("production", "development")
    deployment_options_txt = "deployment can be %r (default is development)" % str.join(', ', deployment_options)
    user_options = [
        ("deployment=", "d", deployment_options_txt),
    ] + InstallCommand.user_options

    def initialize_options(self):
        InstallCommand.initialize_options(self)
        self.deployment = "development"

    def finalize_options(self):
        for deployment_option in self.deployment_options:
            if deployment_option.startswith(self.deployment.lower()):
                self.deployment = deployment_option
                break
        assert self.deployment in self.deployment_options, self.deployment_options_txt
        self.populate_requirements()
        InstallCommand.finalize_options(self)

    def populate_requirements(self):
        reqfn = os.path.join(SourceDir, "requirements.txt")
        requirements = pip.req.parse_requirements(reqfn, session=pip.download.PipSession())
        self.distribution.install_requires = [str(req.req) for req in requirements]

SetupConfiguration = {
    "name": "edge",
    "version": "0.3",
    "author": "Ginkgo Bioworks",
    "author_email": "team@ginkgobioworks.com",
    "description": "Genome Engineering Tool",
    "license": "MIT",
    "packages": ["edge"] + ["edge.%s" % pkg for pkg in setuptools.find_packages("src")],
    "package_dir": {"edge": "src"},
    "include_package_data": True,
    "zip_safe": False,
    "cmdclass": {
        "install": EdgeInstallCommand,
        "install_bower": BowerBuildCommand,
    },
}

def gather_scripts():
    root = os.path.relpath(SourceDir)
    return glob.glob(os.path.join(root, "scripts", "edge_*"))

def gather_packages():
    src_dir = os.path.join(SourceDir, "src")
    packages = setuptools.find_packages(src_dir)
    return ["edge"] + ["edge.%s" % pkg for pkg in packages]

def gather_package_data(config):
    top_level = ["src/templates", "src/static"]
    pkgdata = []
    for dirname in top_level:
        for root, dirs, files in os.walk(dirname):
            root = os.path.relpath(root, "src/")
            pkgdata += [os.path.join(root, fn) for fn in files]
    config["package_data"] = {"edge": pkgdata}

if __name__ == "__main__":
    # allow setup.py to be run from any path
    os.chdir(SourceDir)
    conf = SetupConfiguration.copy()
    conf["scripts"] = gather_scripts()
    conf["packages"] = gather_packages()
    setuptools.setup(**conf)
