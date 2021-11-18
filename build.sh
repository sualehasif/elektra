#!/bin/bash

# check if build exists otherwise make it
if [ ! -d "build" ]; then
    mkdir build
fi

cd build
cmake ..
make
cd ..

