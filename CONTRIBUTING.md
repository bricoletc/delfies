# Reporting bugs and requesting features

If you find a bug with `delfies`, or have a suggested feature, please report them
via the [issue tracker](https://github.com/bricoletc/delfies/issues). Templates 
are provided for bugs and features, but feel free to open a blank issue as well.

# Code contributions: developer instructions

All contributions are warmly welcome!

Before contributing a feature or bugfix, please state this in an open issue
through the issue tracker first, to coordinate with others.

After forking the repository and making changes, open a 
[pull request](https://github.com/bricoletc/delfies/pulls) so we can review it together.

## Checks before opening pull requests

Before creating the pull request, please check that your code runs past
formatting, linting and testing. Formatting relies on
[`black`](https://github.com/psf/black) and
[`isort`](https://github.com/PyCQA/isort), linting relies on
[`flake8`](https://github.com/PyCQA/flake8) and testing relies on
[`pytest`](https://docs.pytest.org/en/stable/). 

To run all these, I've created a Makefile at the root of `delfies`. If you
have the Python packaging tool
[`poetry`](https://github.com/python-poetry/poetry), you can run checks as
follows:

```sh
cd <delfies_directory>
poetry install 
make fmt lint test
# Or 
make precommit
```

## Writing tests

`delfies` has a set of unit and functional/integration tests in the `tests/` directory.
If you add new functions or functionalities, please consider adding corresponding tests.

To check what your new tests cover, you can generate a test coverage report by running:

```sh
cd <delfies_directory>
poetry install 
make coverage
```

This generates a text-file summary and an html report, using [coverage.py](https://github.com/nedbat/coveragepy).
