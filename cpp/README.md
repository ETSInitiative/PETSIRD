# PETSIRD basic C++ example

This directory contains some C++ example code to read/write PETSIRD data. You need to `yardl generate` in the `model` directory first.

The C++ code shows writing to and reading from an HDF5 file

1. Compile the code:

   ```sh
   mkdir -p build
   cmake -G Ninja -S . -B build -DHDF5_ROOT=$CONDA_PREFIX
   ninja
   ```

   If you did not use `conda` to install HDF5, do not add the `-DHDF5_ROOT=$CONDA_PREFIX` part of the `cmake` line.

2. Run the generator: `./petsird_generator test.h5`
3. Run the analyzer: `./petsird_analysis test.h5`
4. You can inspect the HDF5 file by running `h5dump test.h5`
