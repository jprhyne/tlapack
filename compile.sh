#!/bin/env bash
set -x                  # Echo the commands
cmake -B build -GNinja  # Generate the files with Ninja
cmake --build build     # Actually compile
