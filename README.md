# PETSIRD draft data model and examples

[![CI](https://github.com/ETSInitiative/PETSIRD/actions/workflows/ci.yml/badge.svg)](https://github.com/ETSInitiative/PETSIRD/actions/workflows/ci.yml)
[![pre-commit.ci status](https://results.pre-commit.ci/badge/github/ETSInitiative/PETSIRD/main.svg)](https://results.pre-commit.ci/latest/github/ETSInitiative/PETSIRD/main)

The purpose of this repo is to have a working example of a data model for PET imaging (list mode data). This is **not complete**.

## Background

The [Emission Tomography Standardization Initiative (ETSI)](https://etsinitiative.org/)
is working towards establishing a standard for PET Raw Data, called PETSIRD ("PET ETSI Raw Data").

The specification uses the [yardl](https://aka.ms/yardl) tool to define the model.
`yardl` can be used to read the specification (in the `model` directory) and
generate an SDK for both C++ and Python to read/write PETSIRD data.

Browse the [`model`](./model/README.md) directory to get an idea of the content of a PETSIRD file (or stream).

## To get started as a PETSIRD user

If you don't want to modify the yardl model, or modify the PETSIRD code, you probably don't need
this repository at all.

### Python-only users

Just `pip install petsird` or `conda install conda-forge::petsird`. More detail is in the [`python/README.md`](python/README.md).

## C++ users:

You can use `conda install conda-forge::libpetsird` and copy some example code from this repo. More detail is in the [`cpp/README.md`](cpp/README.md).
Alternatively, use https://github.com/ETSInitiative/PETSIRDUseCaseTemplate.

## To get started as a developer

### Quick start

Open this repo in [GitHub Codespaces](https://code.visualstudio.com/docs/remote/codespaces) or
in a [VS Code devcontainer](https://code.visualstudio.com/docs/devcontainers/containers).
This codespace/container will contain all necessary tools, including `yardl` itself, as well as the current repository.

### Do it yourself

Note: these instructions should also work on Windows when you installed Git for Windows with Git Bash and are using a Powershell prompt.

1. create a conda environment with almost all you need and install yardl
   ```sh
   cd /wherever/PETSIRD
   conda env create
   conda activate petsird
   ./install_yardl.sh`
   ```

2. Use `just` to build everything for you
   ```sh
   PATH=$PWD/yardl:$PATH # or in PS $Env:Path = $PWD.Path + "\yardl;" + $Env:Path + $Env:Path
   just build
   ```
   If you cannot use `just` and/or `conda`, you can generate the model yourself.
   ```sh
   cd model
   ../yardl/yardl generate
   cd ..
   ```

3. Have a look at (and try!) the examples in the [`cpp`](cpp/README.md) and/or
[`python`](python/README.md) directories. These README files also contain instructions
for building without `just` and/or `conda`.

Good luck!

## PETSIRD maintainer notes

- Use `pre-commit` before you push to GitHub.
- If it is necessary to update `yardl`, edit the `install_yardl.sh` script.
- When creating a new release, you need to

  1. update the version in `cpp/CMakeLists.txt` (sorry)
  2. merge to `main`
  3. run
  ```sh
  VER=vx.x.x
  git tag -a ${VER} -m ${VER}
  git push origin ${VER}
  ```
  4. edit the draft release automatically created on GitHub and publish it
