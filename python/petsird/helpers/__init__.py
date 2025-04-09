"""
Preliminary helpers for PETSIRD data
"""

#  Copyright (C) 2024 University College London
#
#  SPDX-License-Identifier: Apache-2.0

import typing
from dataclasses import dataclass

import petsird


def get_num_det_els(scanner_geometry: petsird.ScannerGeometry) -> int:
    """Compute total number of detecting elements in the scanner"""
    num_det_els = 0
    for rep_module in scanner_geometry.replicated_modules:
        det_els = rep_module.object.detecting_elements
        for rep_volume in det_els:
            num_det_els += len(rep_volume.transforms) * len(
                rep_module.transforms)
    return num_det_els


@dataclass
class ExpandedDetectionBin:
    """
    Dataclass to store a module ID and element ID, where the latter runs
    over all detecting volumes in the module
    """

    module: int
    el: int
    energy_idx: int


def expand_detection_bins(
    scanner_geometry: petsird.ScannerGeometry,
    detection_bins: typing.Iterable[petsird.DetectionBin]
) -> list[ExpandedDetectionBin]:
    """Find ExpandedDetectionBin for a list of detection_bins"""
    assert len(scanner_geometry.replicated_modules) == 1
    rep_module = scanner_geometry.replicated_modules[0]
    assert len(rep_module.object.detecting_elements) == 1
    num_el_per_module = len(rep_module.object.detecting_elements[0].transforms)

    return [
        ExpandedDetectionBin(module=bin.det_el_idx // num_el_per_module,
                             el=bin.det_el_idx % num_el_per_module,
                             energy_idx=bin.energy_idx)
        for bin in detection_bins
    ]


def get_detection_efficiency(scanner: petsird.ScannerInformation,
                             event: petsird.CoincidenceEvent) -> float:
    """Compute the detection efficiency for a coincidence event"""
    eff = 1

    # per detection_bin efficiencies
    detection_bin_efficiencies = scanner.detection_efficiencies.detection_bin_efficiencies
    if detection_bin_efficiencies is not None:
        eff *= (
            detection_bin_efficiencies[event.detection_bins[0].det_el_idx,
                                       event.detection_bins[0].energy_idx] *
            detection_bin_efficiencies[event.detection_bins[1].det_el_idx,
                                       event.detection_bins[1].energy_idx])

    # per module-pair efficiencies
    module_pair_efficiencies_vector = (
        scanner.detection_efficiencies.module_pair_efficiencies_vector)
    if module_pair_efficiencies_vector is not None:
        module_pair_SGID_LUT = scanner.detection_efficiencies.module_pair_sgidlut
        assert module_pair_SGID_LUT is not None
        expanded_det_bins = expand_detection_bins(scanner.scanner_geometry,
                                                  event.detection_bins)
        assert len(scanner.scanner_geometry.replicated_modules) == 1
        SGID = module_pair_SGID_LUT[expanded_det_bins[0].module,
                                    expanded_det_bins[1].module]
        if SGID < 0:
            return 0.
        module_pair_efficiencies = module_pair_efficiencies_vector[SGID]
        assert module_pair_efficiencies.sgid == SGID
        eff *= module_pair_efficiencies.values[
            expanded_det_bins[0].el,
            event.detection_bins[0].energy_idx,
            expanded_det_bins[1].el,
            event.detection_bins[1].energy_idx,
        ]

    return eff
