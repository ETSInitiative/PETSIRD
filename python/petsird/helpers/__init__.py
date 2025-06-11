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
        num_det_els += len(det_els.transforms) * len(rep_module.transforms)
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
    module_type: petsird.TypeOfModule,
    detection_bins: typing.Iterable[petsird.DetectionBin]
) -> list[ExpandedDetectionBin]:
    """Find ExpandedDetectionBin for a list of detection_bins"""
    rep_module = scanner_geometry.replicated_modules[module_type]
    num_el_per_module = len(rep_module.object.detecting_elements.transforms)

    return [
        ExpandedDetectionBin(module=bin.det_el_idx // num_el_per_module,
                             el=bin.det_el_idx % num_el_per_module,
                             energy_idx=bin.energy_idx)
        for bin in detection_bins
    ]


def expand_detection_bin(
        scanner_geometry: petsird.ScannerGeometry,
        module_type: petsird.TypeOfModule,
        detection_bin: petsird.DetectionBin) -> ExpandedDetectionBin:
    # TODO probably slow implementation, but avoids re-implementation for now
    return expand_detection_bins(scanner_geometry, module_type,
                                 [detection_bin])[0]


def get_detection_efficiency(scanner: petsird.ScannerInformation,
                             module_type_pair: petsird.TypeOfModulePair,
                             event: petsird.CoincidenceEvent) -> float:
    """Compute the detection efficiency for a coincidence event"""
    eff = 1

    # per detection_bin efficiencies
    if scanner.detection_efficiencies is not None:
        detection_bin_efficiencies0 = (
            scanner.detection_efficiencies.detection_bin_efficiencies[
                module_type_pair[0]])
        detection_bin_efficiencies1 = (
            scanner.detection_efficiencies.detection_bin_efficiencies[
                module_type_pair[1]])
        eff *= (
            detection_bin_efficiencies0[event.detection_bins[0].det_el_idx,
                                        event.detection_bins[0].energy_idx] *
            detection_bin_efficiencies1[event.detection_bins[1].det_el_idx,
                                        event.detection_bins[1].energy_idx])

    # per module-pair efficiencies
    module_pair_efficiencies_vectors = (
        scanner.detection_efficiencies.module_pair_efficiencies_vectors)
    if module_pair_efficiencies_vectors is not None:
        module_pair_SGID_LUT = scanner.detection_efficiencies.module_pair_sgidlut
        assert module_pair_SGID_LUT is not None
        expanded_det_bin0 = expand_detection_bin(scanner.scanner_geometry,
                                                 module_type_pair[0],
                                                 event.detection_bins[0])
        expanded_det_bin1 = expand_detection_bin(scanner.scanner_geometry,
                                                 module_type_pair[1],
                                                 event.detection_bins[1])
        assert len(scanner.scanner_geometry.replicated_modules) == 1
        SGID = module_pair_SGID_LUT[module_type_pair[0]][module_type_pair[1]][
            expanded_det_bin0.module, expanded_det_bin1.module]
        if SGID < 0:
            return 0.
        module_pair_efficiencies = module_pair_efficiencies_vectors[
            module_type_pair[0]][module_type_pair[1]][SGID]
        assert module_pair_efficiencies.sgid == SGID
        eff *= module_pair_efficiencies.values[
            expanded_det_bin0.el,
            expanded_det_bin0.energy_idx,
            expanded_det_bin1.el,
            expanded_det_bin1.energy_idx,
        ]

    return eff
