#!/bin/bash

make clean
sphinx-apidoc --force --separate -o apidoc ../xopto
make html
