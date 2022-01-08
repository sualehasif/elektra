#!/bin/bash

# check if build exists otherwise make it
if [ -d "build" ]; then
    rm -rf build
fi
if [ ! -d "build" ]; then
    mkdir build
fi
cd build
cmake -DCMAKE_BUILD_TYPE=Debug .. 
make -j 4
cd ..
