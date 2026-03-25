# PETSIRD basic C++ example using PETSIRD

This directory contains some CMake/C++ example code on how to use
the PETSIRD library as an "external" user.

## Usage:

Assuming you installed PETSIRD in `~/install`, you should be able to
run the following from this directory:
```sh
 cmake -B build  -S . -DCMAKE_PREFIX_PATH=~/install
 cmake  --build build --config Release
 ```

If you have multiple versions of PETSIRD, you can be more specific:
```sh
 cmake -B build -S . -DPETSIRD_DIR=~/install/lib/cmake/PETSIRD-0.7
 ```

