#!/bin/bash

# check if build exists otherwise make it
if [[ -d "build" ]]; then
	rm -rf build
fi
if [[ ! -d "build" ]]; then
	mkdir build
fi
cd build || exit
# cmake -DCMAKE_BUILD_TYPE=Debug -DCILK=1 ..
cmake -DCMAKE_BUILD_TYPE=Debug -DCMAKE_EXPORT_COMPILE_COMMANDS=ON -DSEQUENTIAL=1 ..
make -j 4
cd ..
