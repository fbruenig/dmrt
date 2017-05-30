#!/usr/bin/env bash

#make clean
#make

rm pydmrt.cpython*
python3 setup-pydmrt.py build_ext --inplace
rm build -rf

