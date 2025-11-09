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
    # Post-gen patch: update old xtensor include paths to new layout
    gen="cpp/generated/yardl/detail/ndarray/impl.h"; \
    if [[ -f "$gen" ]]; then \
      perl -0777 -pe \
        's#<xtensor/xarray\.hpp>#<xtensor/containers/xarray.hpp>#g; \
         s#<xtensor/xview\.hpp>#<xtensor/views/xview.hpp>#g; \
         s#<xtensor/xio\.hpp>#<xtensor/io/xio.hpp>#g' \
        -i "$gen"; \
    fi

@build-cpp: generate ensure-configured
    cd cpp/build && cmake --build .

@build-python: generate
    pip install --editable ./python

@build: build-cpp build-python

@run: run-cpp run-python

@run-cpp: build-cpp
    #!/usr/bin/env bash
    cd cpp/build/helpers
    ./petsird_generator testdata.petsird
    ./petsird_analysis testdata.petsird
    rm -f testdata.petsird

@run-python: build-python
    #!/usr/bin/env bash
    python -m petsird.helpers.generator | python -m petsird.helpers.analysis
