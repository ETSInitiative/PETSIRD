# Toy PET data model

The purpose of this repo is to have a simple working example of a data model for PET imaging (list mode data). This is not a working example of how one would actually represent PET data; it is an illustration of how one would use the [yardl](https://github.com/Microsoft/yardl) tool to start defining such a model.

To get started quickly:
1. Open this repo in [GitHub Codespaces](https://code.visualstudio.com/docs/remote/codespaces) or in a [VS Code devcontainer](https://code.visualstudio.com/docs/devcontainers/containers).
2. There is a model package in the `model` directory. `cd` into it and run `yardl generate` to generate C++ code for the model.
3. Compile the code:
    - `cd ../cpp && mkdir -p build && cd build`
    - `cmake ..`
    - `ninja`
3. Run the generator: `./prd_generator`
4. Run the analyzer: `./prd_analysis test.h5`
5. You can inspect the HDF5 file by running `h5dump test.h5`
