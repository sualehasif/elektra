#!/bin/bash

# check if build exists otherwise make it
rm -r build
if [[ ! -d "build" ]]; then
	mkdir build
fi

cd build || exit
cmake -DCMAKE_EXPORT_COMPILE_COMMANDS=ON ..
make -j 4
cd ..
