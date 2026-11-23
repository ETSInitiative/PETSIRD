# PETSIRD C++ library and example

This directory builds a library from C++ code to read/write PETSIRD data.

<strong>Unless you are a PETSIRD developer, you probably want to just install `conda-forge::libppetsird`</strong>.
You can then create your own project based on the [example directory](example/README.md) and the code in the `helpers` directory.

The C++ `helpers` code shows writing to and reading from a PETSIRD file (either in HDF5 or binary format).

These instructions here assume you are familiar with the necessary tools.

> [!TIP]
> If you use conda to install all dependencies, the `justfile` in the main directory will do all of the steps below for you.

## Install dependencies

We support both conda and vcpkg to manage dependencies. However, you will need to get `yardl` yourself (see the [main README](../README.md)).

### Using conda-forge (recommended)

If you don't have conda yet, install miniforge. Then in the main PETSIRD directory, do
```sh
cd PETSIRD
conda env create
conda activate petsird
```
This will give you all tools include CMake, a compatible compiler, as well as all dependencies.

### Using vcpkg (experimental)
Install vcpkg, set `VCPKG_ROOT` as recommended, install CMake, ninja, a compiler and maybe other things, e.g.
on Linux you might need
```sh
sudo apt install pkg-config ninja-build gfortran cmake build-essentials
```
Note: dependencies listed in vcpkg.json will automatically be installed when running CMake as below.
If not, you can always do `vcpkg install`

#### Example configuration on Windows using VS

Install Git for Windows (which includes Git Bash), CMake and Visual Studio
```powershell
winget install Git.Git -e
winget install Kitware.Cmake --override "ADD_CMAKE_TO_PATH=User" --silent
winget install Microsoft.VisualStudio.Community -e --override "--wait --add Microsoft.VisualStudio.Workload.NativeDesktop"
```
For VS, KT selected the "MSVC Build tools", "vcpkg package manager", "C++ CMake Tools", "Windows 11 SDK" and "MSVC v143 - VS2022 C++ build tools".

```powershell
git clone https://github.com/ETSInitiative/PETSIRD
cd PETSIRD
.\install_yardl.sh
cd model
..\yardl\yardl.exe generate
cd ..
```

## Build

You need to `yardl generate` in the `model` directory first.

The instructions below assume your current working directory is `PETSIRD/cpp`.

### Configure the build with CMake

We recommend using one of the presets in the supplied `CMakePresets.json`, e.g.
```sh
cmake --preset conda
```
You can still change values of CMake variables, e.g.
```sh
cmake --preset vcpkg -DCMAKE_INSTALL_PREFIX=/somewhere/nice
```
Alternatively, you can handle things yourself, e.g.
```sh
cmake -S . -B build -DHDF5_ROOT=$CONDA_PREFIX -DCMAKE_INSTALL_PREFIX=~/install
```

You might not need HDF5 support in PETSIRD, in which case you can do
```sh
cmake <other options> -DUSE_HDF5:BOOL=OFF
```
`petsird_generator` will then always output in the PETSIRD binary format.

### Compile the code and install
Adapt the following to your needs
```sh
BUILDDIR=build # or build-vcpkg
cmake --build  ${BUILDDIR} --config Release
cmake --install ${BUILDDIR} --config Release
```

## Try the examples

Assuming you enabled the HDF5 support in PETSIRD during building.

1. Find where the executables are generated. This will be `build/helpers` or `build-vcpkg/helpers/Release` or similar.
2. Run the generator: `<builddir>/petsird_generator test.h5`
3. Run the analyzer: `<builddir>/petsird_analysis test.h5`
4. You can inspect the HDF5 file by running `h5dump test.h5`

You can use the `--binary` option of `petsird_generator` to write in the PETSIRD binary format instead.
`petsird_analysis` will work with both.

## Using PETSIRD as a library

See the [example directory](example/README.md)
