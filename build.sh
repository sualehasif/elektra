#!/bin/bash

# check if build exists otherwise make it
rm -r build
if [ ! -d "build" ]; then
    mkdir build
fi

cd build
cmake -DCILK=1 ..
make -j 4
cd ..

