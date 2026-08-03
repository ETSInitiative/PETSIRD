"""
Preliminary helpers for creating data structures for PETSIRD data
"""

#  Copyright (C) 2026 University College London
#
#  SPDX-License-Identifier: Apache-2.0

from typing import Any

import petsird


def construct_vector(size: int, dtype: Any = float, value=None) -> list:
    """Helper function to create a list of items of type dtype

    Elements will be constructed via the default constructors, so you will still have
    to fill in the actual values.
    """
    if value is None:
        return [dtype() for i in range(size)]
    else:
        return [value for i in range(size)]


def construct_rectangular_matrix(size0: int,
                                 size1: int,
                                 dtype: Any = float,
                                 value: Any = None) -> list[list]:
    """Helper function to create a list of lists of items of type dtype

    The result is equivalent to a 2D array of size (size0, size1).

    Elements will be constructed via the default constructors, so you will still have
    to fill in the actual values.
    """
    return [
        construct_vector(size1, dtype=dtype, value=value) for i in range(size0)
    ]


def construct_lower_triangular_matrix(size: int,
                                      dtype: Any = float,
                                      value: Any = None) -> list[list]:
    """Helper function to create a nested vector equivalent to a lower-triangular matrix

    Result is a 'matrix' where `matrix[i][j]` is only defined/stored for `i>=j`.

    Elements will be constructed via the default constructors, so you will still have
    to fill in the actual values.
    """
    return [
        construct_vector(i + 1, dtype=dtype, value=value) for i in range(size)
    ]


def construct_lower_triangular_or_rectangular_matrix(
        size0: int,
        size1: int,
        is_lower_triangular: bool,
        dtype: Any = float,
        value: Any = None) -> list[list]:
    """Helper function to create a nested vector

    The result is equivelent to either a lower-triangular matrix of size `size0`
    or a matrix size `(size0,size1)`.

    Elements will be constructed via the default constructors, so you will still have
    to fill in the actual values.
    """
    return construct_lower_triangular_matrix(
        size0, dtype=dtype,
        value=value) if is_lower_triangular else construct_rectangular_matrix(
            size0, size1, dtype=dtype, value=value)


def initialize_scanner_information_dimensions(
        scanner: petsird.ScannerInformation, num_module_types: int,
        allocate_detection_bin_efficiencies: bool,
        allocate_module_pair_efficiencies: bool) -> None:
    """Set various structures to have the correct size for the given num_types_of_modules

    This will set scanner.tof_bin_edges, scanner.tof_resolution,
    scanner.event_energy_bin_edges, scanner.energy_resolution_at_511, and
    and (optionally) scanner.detection_efficiencies.detection_bin_efficiencies,
    scanner.detection_efficiencies.module_pair_sgidlut and
    scanner.detection_efficiencies.module_pair_efficiencies_vectors
    to (nested) vectors of the appropriate type and size.

    Elements will be constructed via the default constructors, so you will still have
    to fill in the actual values.

    To prevent dramatic errors, the calibration_factor is set to 1, but this factor has
    to be set correctly afterwards.
    """
    scanner.tof_bin_edges = construct_rectangular_matrix(
        num_module_types, num_module_types, dtype=petsird.BinEdges)
    scanner.tof_resolution = construct_rectangular_matrix(num_module_types,
                                                          num_module_types,
                                                          dtype=float)
    scanner.event_energy_bin_edges = construct_vector(num_module_types,
                                                      dtype=petsird.BinEdges)
    scanner.energy_resolution_at_511 = construct_vector(num_module_types)

    scanner.detection_efficiencies = petsird.DetectionEfficiencies()
    # set this to 1 to avoid the case where the user forgets to set it
    scanner.detection_efficiencies.calibration_factor = 1.

    if allocate_detection_bin_efficiencies:
        scanner.detection_efficiencies.detection_bin_efficiencies = construct_vector(
            num_module_types, dtype=petsird.DetectionBinEfficiencies)

    if allocate_module_pair_efficiencies:
        scanner.detection_efficiencies.module_pair_sgidlut = (
            construct_lower_triangular_matrix(num_module_types,
                                              dtype=petsird.ModulePairSGIDLUT))
        scanner.detection_efficiencies.module_pair_efficiencies_vectors = (
            construct_lower_triangular_matrix(
                num_module_types, dtype=petsird.ModulePairEfficienciesVector))
