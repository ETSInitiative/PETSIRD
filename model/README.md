# PETSIRD draft yardl model

This directory contains several `.yml` files that specify the (draft) PETSIRD data model.
- `Protocol.yml` specifies the sequence of data in the stream
- `*Information.yml` etc define various data-types
- `_package.yml` defines the namespace and where the generated code will be placed.

Use `yardl` to generate C++ and Python code for the model:
  ```sh
  yardl generate
  ```
