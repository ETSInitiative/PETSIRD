# PET data model

The purpose of this repo is to have a simple working example of a data model for PET imaging (list mode data). This is not a working example yet for the actual representation of PET data.

The specification uses the [yardl](https://github.com/Microsoft/yardl) tool to define the model.

## To get started quickly:

1. Open this repo in [GitHub Codespaces](https://code.visualstudio.com/docs/remote/codespaces) or in a [VS Code devcontainer](https://code.visualstudio.com/docs/devcontainers/containers).
2. There is a model package in the `model` directory. `cd` into it and run `yardl generate` to generate C++ and Python code for the model.

### C++

The C++ code shows writing to and reading from an HDF5 file

1. Compile the code:
    - `cd ../cpp && mkdir -p build && cd build`
    - `cmake -G Ninja -S ..` (if you installed HDF5 via `conda`, add `-DHDF5_ROOT=$CONDA_PREFIX`)
    - `ninja`
1. Run the generator: `./prd_generator test.h5`
1. Run the analyzer: `./prd_analysis test.h5`
1. You can inspect the HDF5 file by running `h5dump test.h5`


### Python

The Python code shows piping the compact binary format to standard out and
reading it from standard in.

1. From the repo root `cd python`
1. `python prd_generator.py | python prd_analysis.py`
