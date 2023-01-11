#!/bin/env bash
set -x
cmake -B build -GNinja
cmake --build build
