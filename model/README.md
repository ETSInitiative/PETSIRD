# PETSIRD yardl model

This directory contains several `.yml` files that specify the PETSIRD data model.
- `Protocol.yml` specifies the sequence of data in the stream
- `*Information.yml` etc define various data-types. For instance,
  [ScannerInformation](ScannerInformation.yml) defines the PET scanner in terms of
  [detectors](DetectorInformation.yml), [materials](MaterialInformation.yml) and various
  properties such as energy windows etc.
- `_package.yml` defines the namespace and where the generated code will be placed.

Note that `yardl` reads all `.yml` files, so the order of the type definitions
does not matter as far as `yardl` concerns.
You can read more about the yardl syntax here: https://aka.ms/yardl/docs

Use `yardl` to generate C++ and Python code for the model:
  ```sh
  yardl generate
  ```
