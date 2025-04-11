set shell := ['bash', '-ceuo', 'pipefail']

@default: build

@configure:
    mkdir -p cpp/build; \
    cd cpp/build; \
    cmake -GNinja ..

@ensure-configured:
    if [ ! -f cpp/build/CMakeCache.txt ]; then \
        just configure; \
    fi

@generate:
    cd model && yardl generate

@build-cpp: generate ensure-configured
    cd cpp/build && ninja

@build-python: generate
    pip install --editable ./python

@build: build-cpp build-python

@run: run-cpp run-python

@run-cpp: build-cpp
    #!/usr/bin/env bash
    cd cpp/build
    ./petsird_generator testdata.petsird
    ./petsird_analysis testdata.petsird
    rm -f testdata.petsird

@run-python: build-python
    #!/usr/bin/env bash
    python -m petsird.helpers.generator | python -m petsird.helpers.analysis
