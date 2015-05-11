#!/bin/bash
./build.sh
mpirun -np ${1} ./debug
