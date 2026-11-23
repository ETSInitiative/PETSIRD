# PETSIRD C++ library and example

This directory builds a library from C++ code to read/write PETSIRD data. You need to `yardl generate` in the `model` directory first.

The C++ code shows writing to and reading from a PETSIRD file (either in HDF5 or binary format).

1. Compile the code:

   ```sh
   cmake -S cpp -B build -DHDF5_ROOT=$CONDA_PREFIX -DCMAKE_INSTALL_PREFIX=install
   cmake --install build
   ```

   If you did not use `conda` to install HDF5, do not add the `-DHDF5_ROOT=$CONDA_PREFIX` part of the `cmake` line.

2. Run the generator: `install/petsird_generator test.h5`
3. Run the analyzer: `install/petsird_analysis test.h5`
4. You can inspect the HDF5 file by running `h5dump test.h5`

You can use the `--binary` option of `petsird_generator` to write in the PETSIRD binary format instead.
`petsird_analysis` will work with both.

## Disabling HDF5 support

You might not need HDF5 support in PETSIRD, in which case you can do
```sh
cmake <other options> -DUSE_HDF5:BOOL=OFF
```
`petsird_generator` will then always output in the PETSIRD binary format.

## Using this is a library

See the [example directory](example/README.md)
