#!/bin/bash

rm -rf proteus/xtensor
mkdir -p proteus/xtensor
cd proteus/xtensor
mkdir xtl xtensor-python xtensor pybind11
wget https://github.com/QuantStack/xtl/archive/0.6.11.tar.gz            -O xtl.tar.gz
wget https://github.com/QuantStack/xtensor-python/archive/0.24.1.tar.gz -O xtensor-python.tar.gz
wget https://github.com/xtensor-stack/xtensor/archive/0.21.2.tar.gz     -O xtensor.tar.gz
wget https://github.com/pybind/pybind11/archive/v2.4.3.tar.gz           -O pybind11.tar.gz
tar zxf xtl.tar.gz              --strip-components=1 -C xtl
tar zxf xtensor-python.tar.gz   --strip-components=1 -C xtensor-python
tar zxf xtensor.tar.gz          --strip-components=1 -C xtensor
tar zxf pybind11.tar.gz         --strip-components=1 -C pybind11
rm *.gz
