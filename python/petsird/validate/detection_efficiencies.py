"""
Validation functions and helpers related to scanner detection efficiencies.
"""
#  Copyright (C) 2026 University College London
#
#  SPDX-License-Identifier: Apache-2.0

import petsird
from petsird.helpers import indices_in_2D_LUT
from petsird.helpers.geometry import are_possibly_symmetric_with_first


def validate_module_pair_LUTs(scanner: petsird.ScannerInformation) -> None:
    """check SGID values and sizes of LUTs

    Uses `are_possibly_symmetric_with_first` to check if module-pairs with
    the same SGID have the same distances between detecting elements.
    """
    num_module_types = scanner.scanner_geometry.number_of_module_types()
    if scanner.detection_efficiencies.module_pair_sgidlut:
        for mtype0 in range(num_module_types):
            for mtype1 in range(mtype0 + 1):
                LUT = scanner.detection_efficiencies.module_pair_sgidlut[
                    mtype0][mtype1]
                num_SGIDs = max(max(row) for row in LUT) + 1
                if num_SGIDs != len(
                        scanner.detection_efficiencies.
                        module_pair_efficiencies_vectors[mtype0][mtype1]):
                    raise ValueError(
                        "length of module_pair_efficiencies_vectors not equal to max SGID"
                    )

                for SGID in range(num_SGIDs):
                    if not all(
                            are_possibly_symmetric_with_first(
                                scanner, mtype0, mtype1,
                                indices_in_2D_LUT(SGID, LUT))):
                        raise ValueError(
                            f"module_pair_SGIDLUT[{mtype0}][{mtype1}] SGID={SGID} "
                            + "contains pairs with different distances")
