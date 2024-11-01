# PETSIRD draft data model and examples

The purpose of this repo is to have a working example of a data model for PET imaging (list mode data). This is **not complete**, and aspects do not follow the actual proposal of PETSIRD yet.

## Background

The [Emission Tomography Standardization Initiative (ETSI)](https://etsinitiative.org/)
is working towards establishing a standard for PET Raw Data, called PETSIRD ("PET ETSI Raw Data").

The specification uses the [yardl](https://aka.ms/yardl) tool to define the model.
`yardl` can be used to read the specification (in the `model` directory) and
generate an SDK for both C++ and Python to read/write PETSIRD data.

## To get started as a Python user:

If you don't want to modify the yardl model, clone this repository and follow instructions in the [`python/README.md`](python/README.md).

## To get started quickly as a developer:

1. Open this repo in [GitHub Codespaces](https://code.visualstudio.com/docs/remote/codespaces) or
in a [VS Code devcontainer](https://code.visualstudio.com/docs/devcontainers/containers).
This codespace/container will contain all necessary tools, including `yardl` itself, as well as the current repository.
2. Browse the [`model`](./model/README.md) directory.
3. Use `yardl` to generate C++ and Python code for the model:

   ```sh
   cd model
   yardl generate
   ```

4. Have a look at (and try!) the examples in the [`cpp`](cpp/README.md) and/or
[`python`](python/README.md) directories.


