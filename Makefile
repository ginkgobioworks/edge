.PHONY: help all image \
	clean clean-pyc clean-build clean-js \
	${SETUP_COMMANDS} \
	${MANAGE_COMMANDS} watch \
	test-all test-ci \
	bump/major bump/minor bump/patch \
	start \
	release

help:
	@echo "clean-build - remove build artifacts"
	@echo "clean-pyc - remove Python file artifacts"
	@echo "clean - remove both build and Python artifacts"
	@echo "install_bower - install JS dependencies with bower"
	@echo "build_assets - build the static assets"
	@echo "migrate, syncdb - set up the db"
	@echo "add-s288c - load example data into the db"
	@echo "dbshell, shell - run a sql or python shell with the app & db loaded"
	@echo "test - run base python tests"
	@echo "flake8 - check style with flake8"
	@echo "run - run the dev server"
	@echo "start - run migration and run the dev server"
	@echo "watch - watch static assets and recompile automatically"
	@echo "bash - start a bash shell; only really useful when run as bash-ext or bash-ext_fast"
	@echo "sdist, bdist_wheel - package"
	@echo "bump/major bump/minor bump/patch - bump the version"
	@echo "release - package and upload a release"
	@echo "image - make the image"
	@echo "build_edge_blastdb, build_genome_blastdb - custom manage.py django commands"
	@echo "export_gff, import_gff - custom manage.py django commands"
	@echo "remove_fragment, remove_genome - custom manage.py django commands"
	@echo "push - push the image to the docker registry"
	@echo "NB: add '-ext' to any target to run it inside the image from the host machine; will"
	@echo "    rebuild the image."
	@echo "    add '-ext_fast' to any target to run it from the host machine in the image"
	@echo "    without rebuilding"


PROJECT_NAME = edge
BUILD_IMAGE ?= ${PROJECT_NAME}
EDGE_HOME ?= /usr/src/${PROJECT_NAME}

MANAGE = cd src; python manage.py
SETUP = python setup.py

all: test-all


# Python packaging commands from setup.py

SETUP_COMMANDS = install_bower build_assets sdist bdist_wheel flake8
${SETUP_COMMANDS}:
	${SETUP} $@ ${args}


# Cleanup

clean: clean-build clean-pyc clean-js

clean-build:
	rm -rf build/
	rm -rf dist/
	rm -rf src/*.egg-info

clean-pyc:
	find . -name '*.pyc' -exec rm -f {} +
	find . -name '*.pyo' -exec rm -f {} +
	find . -name '*~' -exec rm -f {} +

clean-js:
	rm -rf src/edge/static/edge/lib/*
	rm -vf src/edge/static/edge/edge{.css,.html,.js,.min.js,_jst.js}


# Testing

test-all: test flake8

test-ci: migrate test

# Django management commands from manage.py

MANAGE_COMMANDS = \
  build_edge_blastdb build_genome_blastdb \
  export_gff import_gff \
  remove_fragment remove_genome \
  createsuperuser dbshell shell test

${MANAGE_COMMANDS}:
	${MANAGE} $@ ${args}

start: clean migrate run

run: install_bower build_assets
	${MANAGE} runserver ${EDGE_IP}:${EDGE_PORT}

migrate: syncdb
	${MANAGE} migrate

syncdb:
	sleep 5
	${MANAGE} syncdb --noinput

watch: install_bower
	${MANAGE} assets watch


# Bump the version and release

bump/major bump/minor bump/patch:
	bumpversion --verbose $(@F)

release: clean sdist bdist_wheel
	ls -l dist
	twine upload dist/*


# Utilities

bash:
	bash

# Example data

add-s288c:
	${MANAGE} import_gff 'Saccharomyces cerevisiae' ../example/sc_s288c.gff
	${MANAGE} build_edge_blastdb


# Generically execute make targets from outside the Docker container

MAKE_EXT = docker-compose run --rm --service-ports ${PROJECT_NAME} make -C ${EDGE_HOME}

%-ext : image
	${MAKE_EXT} $*

# Run a command from outside the image without rebuilding the image; useful when you already have
# a version of the image running
%-ext_fast :
	${MAKE_EXT} $*

# Build the image

image:
	GIT_USER_NAME=`git config user.name` GIT_USER_EMAIL=`git config user.email` docker-compose build --pull ${PROJECT_NAME}

#  vim: set noet :
