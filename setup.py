#!/usr/bin/env python
from __future__ import print_function

import os
from subprocess import check_output, STDOUT

from setuptools import setup, find_packages, Command
from setuptools.command.sdist import sdist


class BowerInstallCommand(Command):
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
    print(check_output(cmd, stderr=STDOUT))


class BuildAssetsCommand(Command):
  description = 'build django assets'
  user_options = []

  def initialize_options(self):
    pass

  def finalize_options(self):
    pass

  def run(self):
    src_dir = os.path.join(os.path.dirname(__file__), 'src')
    print(check_output(['python', 'manage.py', 'assets', 'build'], cwd=src_dir, stderr=STDOUT))


class SdistCommandWithJS(sdist):
  def run(self):
    self.run_command('install_bower')
    self.run_command('build_assets')
    sdist.run(self)


setup(
  name='edge-genome',
  version='0.3.0',

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

  setup_requires=[
    'django_assets',
  ],
  install_requires=[
    'django < 1.7',
    'jsmin',
    'celery < 4.0',
    'django-celery',
    'bcbio-gff',
    'pytz >= 1023b',
  ],
  tests_require=[
    'nose'
    'django-nose < 1.4',
    'flake8',
    'mysql-python',
  ],
)
