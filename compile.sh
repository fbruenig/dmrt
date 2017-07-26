#!/usr/bin/env bash

#make clean
#make

rm __pycache__* -rf
rm difftools.pyc -f
rm pydmrt.so
rm pydmrt.cpython*
python3 setup-pydmrt.py build_ext --inplace
rm build -rf

