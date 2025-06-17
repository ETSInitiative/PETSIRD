"""
Preliminary helpers for PETSIRD data
"""

#  Copyright (C) 2024 - 2025 University College London
#
#  SPDX-License-Identifier: Apache-2.0

import typing
from dataclasses import dataclass

import petsird


def get_num_det_els(scanner_geometry: petsird.ScannerGeometry,
                    type_of_module: petsird.TypeOfModule) -> int:
    """Compute total number of detecting elements in all modules of the given type"""
    rep_module = scanner_geometry.replicated_modules[type_of_module]
    det_els = rep_module.object.detecting_elements
    return len(det_els.transforms) * len(rep_module.transforms)


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
    type_of_module: petsird.TypeOfModule,
    detection_bins: typing.Iterable[petsird.DetectionBin]
) -> list[ExpandedDetectionBin]:
    """Find ExpandedDetectionBin for a list of detection_bins"""
    rep_module = scanner_geometry.replicated_modules[type_of_module]
    num_el_per_module = len(rep_module.object.detecting_elements.transforms)

    return [
        ExpandedDetectionBin(module=bin.det_el_idx // num_el_per_module,
                             el=bin.det_el_idx % num_el_per_module,
                             energy_idx=bin.energy_idx)
        for bin in detection_bins
    ]


def expand_detection_bin(
        scanner_geometry: petsird.ScannerGeometry,
        type_of_module: petsird.TypeOfModule,
        detection_bin: petsird.DetectionBin) -> ExpandedDetectionBin:
    # TODO probably slow implementation, but avoids re-implementation for now
    return expand_detection_bins(scanner_geometry, type_of_module,
                                 [detection_bin])[0]


def get_detection_efficiency(scanner: petsird.ScannerInformation,
                             type_of_module_pair: petsird.TypeOfModulePair,
                             event: petsird.CoincidenceEvent) -> float:
    """Compute the detection efficiency for a coincidence event"""
    eff = 1

    # per detection_bin efficiencies
    if scanner.detection_efficiencies is not None:
        detection_bin_efficiencies0 = (
            scanner.detection_efficiencies.detection_bin_efficiencies[
                type_of_module_pair[0]])
        detection_bin_efficiencies1 = (
            scanner.detection_efficiencies.detection_bin_efficiencies[
                type_of_module_pair[1]])
        eff *= (
            detection_bin_efficiencies0[event.detection_bins[0].det_el_idx,
                                        event.detection_bins[0].energy_idx] *
            detection_bin_efficiencies1[event.detection_bins[1].det_el_idx,
                                        event.detection_bins[1].energy_idx])
        if eff == 0:
            return 0.

    # per module-pair efficiencies
    module_pair_efficiencies_vectors = (
        scanner.detection_efficiencies.module_pair_efficiencies_vectors)
    if module_pair_efficiencies_vectors is not None:
        module_pair_SGID_LUT = scanner.detection_efficiencies.module_pair_sgidlut
        assert module_pair_SGID_LUT is not None
        expanded_det_bin0 = expand_detection_bin(scanner.scanner_geometry,
                                                 type_of_module_pair[0],
                                                 event.detection_bins[0])
        expanded_det_bin1 = expand_detection_bin(scanner.scanner_geometry,
                                                 type_of_module_pair[1],
                                                 event.detection_bins[1])
        assert len(scanner.scanner_geometry.replicated_modules) == 1
        SGID = module_pair_SGID_LUT[type_of_module_pair[0]][
            type_of_module_pair[1]][expanded_det_bin0.module,
                                    expanded_det_bin1.module]
        if SGID < 0:
            return 0.
        module_pair_efficiencies = module_pair_efficiencies_vectors[
            type_of_module_pair[0]][type_of_module_pair[1]][SGID]
        assert module_pair_efficiencies.sgid == SGID
        eff *= module_pair_efficiencies.values[
            expanded_det_bin0.el,
            expanded_det_bin0.energy_idx,
            expanded_det_bin1.el,
            expanded_det_bin1.energy_idx,
        ]

    return eff
