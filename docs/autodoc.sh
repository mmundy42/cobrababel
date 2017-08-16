#!/usr/bin/env bash
rm cobrababel*.rst
sphinx-apidoc -o . ../cobrababel ../cobrababel/test
rm modules.rst
jupyter nbconvert --to=rst *.ipynb