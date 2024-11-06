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

@build: generate ensure-configured
    cd cpp/build && ninja

@run: run-cpp run-python

@run-cpp: build
    #!/usr/bin/env bash
    cd cpp/build
    ./petsird_generator testdata.petsird
    ./petsird_analysis testdata.petsird
    rm -f testdata.petsird

@run-python: generate
    #!/usr/bin/env bash
    cd python
    python petsird_generator.py | python petsird_analysis.py
