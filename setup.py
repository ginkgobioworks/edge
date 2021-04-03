#!/usr/bin/env python
from __future__ import print_function

import os
import sys
from subprocess import (
    CalledProcessError,
    STDOUT,
    check_output,
)

from setuptools import setup, find_packages, Command
from setuptools.command.sdist import sdist


class SubproccessCommand(Command):
    def _run_with_output(self, *args, **kwargs):
        try:
            print(check_output(*args, **kwargs))
        except CalledProcessError as e:
            print(e.output, file=sys.stderr)

            raise e


class BowerInstallCommand(SubproccessCommand):
    description = 'run bower commands [install by default]'
    user_options = [
        ('bower-command=', 'c', 'Bower command to run. Defaults to "install".'),
        ('force-latest', 'F', 'Force latest version on conflict.'),
        ('allow-root', 'R', 'Allow bower to be run as root.'),
        ('production', 'p', 'Do not install project devDependencies.'),
    ]
    boolean_options = ['production', 'force-latest', 'allow-root']

    def initialize_options(self):
        self.force_latest = False
        self.production = False
        self.bower_command = 'install'
        self.allow_root = False

    def finalize_options(self):
        pass

    def run(self):
        cmd = ['bower', self.bower_command, '--no-color']
        if self.force_latest:
            cmd.append('-F')
        if self.production:
            cmd.append('-p')
        if self.allow_root:
            cmd.append('--allow-root')
        self._run_with_output(cmd, stderr=STDOUT)


class BuildAssetsCommand(SubproccessCommand):
    description = 'build django assets'
    user_options = []

    def initialize_options(self):
        pass

    def finalize_options(self):
        pass

    def run(self):
        src_dir = os.path.join(os.path.dirname(__file__), 'src')
        self._run_with_output(['python', 'manage.py', 'assets', 'build'],
                              cwd=src_dir, stderr=STDOUT)


class SdistCommandWithJS(sdist):
    def run(self):
        self.run_command('install_bower')
        self.run_command('build_assets')
        sdist.run(self)


setup(
    name='edge-genome',
    version='2.17.0',

    author='Ginkgo Bioworks',
    author_email='devs@ginkgobioworks.com',
    url='https://github.com/ginkgobioworks/edge/',

    description='Genome Engineering Tool',
    long_description=open('README.rst').read(),

    license='MIT',
    keywords='edge genome annotate track query django database ginkgo',
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Environment :: Web Environment',
        'Environment :: Other Environment',
        'Framework :: Django :: 1.6',
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Operating System :: OS Independent',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: JavaScript',
        'Programming Language :: SQL',
        'Topic :: Database',
        'Topic :: Internet :: WWW/HTTP :: Indexing/Search',
        'Topic :: Internet :: WWW/HTTP :: WSGI :: Application',
        'Topic :: Scientific/Engineering :: Artificial Intelligence',
        'Topic :: Scientific/Engineering :: Artificial Life',
        'Topic :: Scientific/Engineering :: Information Analysis',
        'Topic :: Software Development :: Libraries',
        'Topic :: Software Development :: Libraries :: Python Modules',
        'Topic :: Text Processing :: General',
    ],

    cmdclass={
        'install_bower': BowerInstallCommand,
        'build_assets': BuildAssetsCommand,
        'sdist': SdistCommandWithJS,
    },

    package_dir={'': 'src'},
    packages=find_packages('src', exclude=['*.migrations', '*server*']),
    include_package_data=True,
    zip_safe=True,

    # Django has to be specified here because django_assets tries to install Django>=1.7 and will
    # force install Django2+ in setup
    setup_requires=[
        'django_assets >= 0.12',
        'Django ~= 2.2.6',
    ],
    install_requires=[
        'django ~= 2.2.6',
        'jsmin',
        'celery >= 4.0',
        'bcbio-gff',
        'pytz >= 1023b',
    ],
    tests_require=[
        'nose',
        'coverage',
        'django-nose',
        'flake8',
        'mysql-python',
    ],
)
