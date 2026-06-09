set shell := ['bash', '-ceuo', 'pipefail']

cmake_install_prefix := "$CONDA_PREFIX"
cmake_build_type := "Release"
cmake_build_dir := "cpp/build"

@default: build

@help:
    echo 'Usage: "just [cmake_install_prefix=<CMAKE_INSTALL_PREFIX>] [cmake_build_dir=<BUILD_DIR>] [cmake_build_type=Release]"'
    echo 'The variables will be passed to CMake for the C++ build.'
    echo 'Defaults:'
    echo "  cmake_install_prefix=\$CONDA_PREFIX (which is currently set to \"$CONDA_PREFIX\")"
    echo '  cmake_build_type=Release'
    echo '  cmake_build_dir=cpp/build (this is relative to the PETSIRD folder)'
    echo 'Run "just --summary" for possible recipes (default recipe is "build")'

@configure: generate
    cmake -GNinja -S cpp -B {{cmake_build_dir}} \
      -DCMAKE_BUILD_TYPE:BOOL={{cmake_build_type}} \
      -DCMAKE_INSTALL_PREFIX:PATH={{cmake_install_prefix}}

@ensure-configured:
    # Need to repeat variables here, as we do another call to "just"
    if [ ! -f {{cmake_build_dir}}/CMakeCache.txt ]; then \
      just cmake_install_prefix={{cmake_install_prefix}} \
        cmake_build_type={{cmake_build_type}} cmake_build_dir={{cmake_build_dir}} \
        configure; \
    fi

@generate:
    cd model && yardl generate

@build-cpp: generate ensure-configured
    cd {{cmake_build_dir}} && \
    cmake --build . --config {{cmake_build_type}} && \
    cmake --install .

@build-python: generate
    python -m pip install --editable ./python

@build: build-cpp build-python

@run: run-cpp run-python

@run-cpp: build-cpp
    #!/usr/bin/env bash
    cd {{cmake_build_dir}}/helpers
    ./petsird_generator testdata.petsird
    ./petsird_analysis testdata.petsird
    rm -f testdata.petsird

@run-python: build-python
    #!/usr/bin/env bash
    python -m petsird.helpers.generator | python -m petsird.helpers.analysis
