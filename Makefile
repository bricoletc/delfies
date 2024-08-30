# Strongly inspired from: https://github.com/snakemake/snakefmt
PROJECT = delfies
COVG_REPORT = htmlcov/index.html
OS := $(shell uname -s)
VERSION := $(shell poetry version -s)
BOLD := $(shell tput bold)
NORMAL := $(shell tput sgr0)
# MAIN #########################################################################

.PHONY: all
all: install

# DEPENDENCIES #################################################################
.PHONY: install
install:
	poetry install

.PHONY: install-ci
install-ci:
	poetry install --no-interaction
	poetry run delfies --version

# TIDY #################################################################
.PHONY: fmt
fmt:
	poetry run isort . --skip tmp_work --skip venv --skip dist
	poetry run black .

.PHONY: lint
lint:
	poetry run flake8 . --exclude tmp_work,venv,dist

.PHONY: check-fmt
check-fmt:
	poetry run isort --check-only .
	poetry run black --check .

# TEST ########################################################################
.PHONY: test
test:
	poetry run pytest tests/

.PHONY: coverage
coverage:
	poetry run pytest --cov-report term --cov-report html --cov=$(PROJECT) --cov-branch
ifeq ($(OS), Linux)
	xdg-open $(COVG_REPORT)
else ifeq ($(OS), Darwin)
	open $(COVG_REPORT)
else
	echo "ERROR: Unknown OS detected - $OS"
endif

# PRECOMMIT ########################################################################
# runs format, lint and test
.PHONY: precommit
precommit: fmt lint test

# BUILD ########################################################################
.PHONY: build
build:
	poetry build

# TAG ########################################################################
# prints out the commands to run to tag the release and push it
.PHONY: tag
tag:
	@echo "Run $(BOLD)git tag -a $(VERSION) -m <message>$(NORMAL) to tag the release"
	@echo "Then run $(BOLD)git push upstream $(VERSION)$(NORMAL) to push the tag"
