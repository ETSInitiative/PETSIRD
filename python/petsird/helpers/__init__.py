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
class ModuleAndElement:
    """
    Dataclass to store a module ID and element ID, where the latter runs
    over all detecting volumes in the module
    """

    module: int
    el: int


def get_module_and_element(
        scanner_geometry: petsird.ScannerGeometry,
        scanner_det_ids: typing.Iterable[int]) -> list[ModuleAndElement]:
    """Find ModuleAndElement for a list of detector_ids"""
    assert len(scanner_geometry.replicated_modules) == 1
    rep_module = scanner_geometry.replicated_modules[0]
    assert len(rep_module.object.detecting_elements) == 1
    num_el_per_module = len(rep_module.object.detecting_elements[0].ids)

    return [
        ModuleAndElement(module=det // num_el_per_module,
                         el=det % num_el_per_module) for det in scanner_det_ids
    ]


def get_detection_efficiency(scanner: petsird.ScannerInformation,
                             event: petsird.CoincidenceEvent) -> float:
    """Compute the detection efficiency for a coincidence event"""
    eff = 1

    # per det_el efficiencies
    det_el_efficiencies = scanner.detection_efficiencies.det_el_efficiencies
    if det_el_efficiencies is not None:
        eff *= (det_el_efficiencies[event.detector_ids[0],
                                    event.energy_indices[0]] *
                det_el_efficiencies[event.detector_ids[1],
                                    event.energy_indices[1]])

    # per module-pair efficiencies
    module_pair_efficiencies_vector = (
        scanner.detection_efficiencies.module_pair_efficiencies_vector)
    if module_pair_efficiencies_vector is not None:
        module_pair_SGID_LUT = scanner.detection_efficiencies.module_pair_sgidlut
        assert module_pair_SGID_LUT is not None
        mod_and_els = get_module_and_element(scanner.scanner_geometry,
                                             event.detector_ids)
        assert len(scanner.scanner_geometry.replicated_modules) == 1
        SGID = module_pair_SGID_LUT[mod_and_els[0].module,
                                    mod_and_els[1].module]
        assert SGID >= 0
        module_pair_efficiencies = module_pair_efficiencies_vector[SGID]
        assert module_pair_efficiencies.sgid == SGID
        eff *= module_pair_efficiencies.values[
            mod_and_els[0].el,
            event.energy_indices[0],
            mod_and_els[1].el,
            event.energy_indices[1],
        ]

    return eff
