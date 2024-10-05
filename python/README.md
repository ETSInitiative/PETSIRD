# PETSIRD basic Python example

This directory contains some Python example code to read/write PETSIRD data. You need to `yardl generate` in the `model` directory first.

The Python code shows piping the compact binary format to standard out and
reading it from standard in. This can be used as follows:

1. From the repo root `cd python`
1. `python prd_generator.py | python prd_analysis.py`

There is also a very basic utility to plot the scanner geometry. For instance
```
python prd_generator.py > test.bin
python prd_plot_scanner.py < test.bin`

