#!/bin/env bash
set -x                  # Echo the commands
cmake -B build -GNinja  # Generate the files with Ninja. 
#To enable debugging symbols replace above line with cmake -DCMAKE_BUILD_TYPE=DEBUG --build build
cmake --build build     # Actually compile
#/home/langou/Repositories/tlapack.git/build/test/test_eispack_hqr
