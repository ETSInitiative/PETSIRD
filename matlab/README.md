# PETSIRD basic MATLAB example

This directory contains some MATLAB example code to read/write PETSIRD data.

> WARNING: This code is likely out-of-date with the rest of this repository.
> It is untested. Hopefully, it will occasionally be updated.

## Content

- Generator and Analysis tools ported from Python to MATLAB (`generator.m` and `analysis.m`)
- `matlab/build-toolbox.sh` script will build a portable MATLAB toolbox containing the tools

## Testing

From the top-level of the repository, the following script exercises the MATLAB tools, the output of which can be visually compared to the output of the Python tools.
```bash
#!/usr/bin/env bash

matlab -batch "addpath('matlab/toolbox/'); petsird.helpers.generator('testdata.matlab.petsird');"
python  -m petsird.helpers.generator.py > testdata.python.petsird

matlab -batch "addpath('matlab/toolbox/'); petsird.helpers.analysis('testdata.matlab.petsird');" > matlab2matlab.analysis
matlab -batch "addpath('matlab/toolbox/'); petsird.helpers.analysis('testdata.python.petsird');" > python2matlab.analysis

cat testdata.python.petsird | python -m petsird.helpers.analysis.py > python2python.analysis
cat testdata.matlab.petsird | python -m petsird.helpers.analysis.py > matlab2python.analysis
```
