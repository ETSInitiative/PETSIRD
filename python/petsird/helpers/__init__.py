"""
Preliminary helpers for PETSIRD data
"""

#  Copyright (C) 2024 - 2025 University College London
#
#  SPDX-License-Identifier: Apache-2.0

import typing

import petsird


# TODO remove?
def get_num_det_els(scanner: petsird.ScannerInformation,
                    type_of_module: petsird.TypeOfModule) -> int:
    """Compute total number of detecting elements in all modules of the given type"""
    rep_module = scanner.scanner_geometry.replicated_modules[type_of_module]
    det_els = rep_module.object.detecting_elements
    return len(det_els.transforms) * len(rep_module.transforms)


def get_num_detection_bins(scanner: petsird.ScannerInformation,
                           type_of_module: petsird.TypeOfModule) -> int:
    """Compute total number of detecting bins in all modules of the given type"""
    rep_module = scanner.scanner_geometry.replicated_modules[type_of_module]
    energy_bin_edges = scanner.event_energy_bin_edges[type_of_module]
    det_els = rep_module.object.detecting_elements
    return len(det_els.transforms) * len(
        rep_module.transforms) * energy_bin_edges.number_of_bins()


def expand_detection_bins(
    scanner: petsird.ScannerInformation, type_of_module: petsird.TypeOfModule,
    detection_bins: typing.Iterable[petsird.DetectionBin]
) -> list[petsird.ExpandedDetectionBin]:
    """Find ExpandedDetectionBin for a list of detection_bins"""
    rep_module = scanner.scanner_geometry.replicated_modules[type_of_module]
    num_el_per_module = len(rep_module.object.detecting_elements.transforms)
    energy_bin_edges = scanner.event_energy_bin_edges[type_of_module]
    num_en = energy_bin_edges.number_of_bins()

    return [
        petsird.ExpandedDetectionBin(
            module_index=bin // (num_el_per_module * num_en),
            element_index=(bin // num_en) % num_el_per_module,
            energy_index=bin % num_en) for bin in detection_bins
    ]


def expand_detection_bin(
        scanner: petsird.ScannerInformation,
        type_of_module: petsird.TypeOfModule,
        detection_bin: petsird.DetectionBin) -> petsird.ExpandedDetectionBin:
    # TODO probably slow implementation, but avoids re-implementation for now
    return expand_detection_bins(scanner, type_of_module, [detection_bin])[0]


def make_detection_bins(
    scanner: petsird.ScannerInformation, type_of_module: petsird.TypeOfModule,
    expanded_detection_bins: typing.Iterable[petsird.ExpandedDetectionBin]
) -> list[petsird.DetectionBin]:
    """Find DetectionBin for a list of expanded_detection_bins"""
    rep_module = scanner.scanner_geometry.replicated_modules[type_of_module]
    num_el_per_module = len(rep_module.object.detecting_elements.transforms)
    energy_bin_edges = scanner.event_energy_bin_edges[type_of_module]
    num_en = energy_bin_edges.number_of_bins()

    return [
        (bin.energy_index +
         (bin.element_index + bin.module_index * num_el_per_module) * num_en)
        for bin in expanded_detection_bins
    ]


def make_detection_bin(
    scanner: petsird.ScannerInformation, type_of_module: petsird.TypeOfModule,
    expanded_detection_bin: petsird.ExpandedDetectionBin
) -> petsird.DetectionBin:
    """Find DetectionBin for an expanded_detection_bin"""
    # TODO probably slow implementation, but avoids re-implementation for now
    return make_detection_bins(scanner, type_of_module,
                               [expanded_detection_bin])[0]


@typing.overload
def get_detection_efficiency(scanner: petsird.ScannerInformation,
                             type_of_module_pair: petsird.TypeOfModulePair,
                             detection_bin_1: petsird.DetectionBin,
                             detection_bin_2: petsird.DetectionBin,
                             with_calibration_factor: bool = ...) -> float:
    """Compute the detection efficiency for a pair of detectors"""
    ...


@typing.overload
def get_detection_efficiency(scanner: petsird.ScannerInformation,
                             type_of_module_pair: petsird.TypeOfModulePair,
                             event: petsird.CoincidenceEvent,
                             with_calibration_factor: bool = ...) -> float:
    """Compute the detection efficiency for a coincidence event"""
    ...


_DetectionBinUnannotated = typing.get_args(petsird.DetectionBin)[0]


def get_detection_efficiency(scanner: petsird.ScannerInformation,
                             type_of_module_pair: petsird.TypeOfModulePair,
                             event_or_detection_bin_1: typing.Union[
                                 petsird.CoincidenceEvent,
                                 petsird.DetectionBin],
                             detection_bin_2: petsird.DetectionBin = None,
                             with_calibration_factor: bool = True) -> float:
    """Compute the detection efficiency"""
    if isinstance(event_or_detection_bin_1, _DetectionBinUnannotated):
        detection_bin_1 = event_or_detection_bin_1
        assert detection_bin_2 is not None, "Second detection bin must be provided"
    else:
        detection_bin_1, detection_bin_2 = event_or_detection_bin_1.detection_bins[:
                                                                                   2]

    if scanner.detection_efficiencies is None:
        # should never happen really, but this way, we don't crash.
        return 1.

    eff = (scanner.detection_efficiencies.calibration_factor
           if with_calibration_factor else 1.)

    # per detection_bin efficiencies
    detection_bin_efficiencies = scanner.detection_efficiencies.detection_bin_efficiencies
    if detection_bin_efficiencies is not None:
        detection_bin_efficiencies0 = (
            detection_bin_efficiencies[type_of_module_pair[0]])
        detection_bin_efficiencies1 = (
            detection_bin_efficiencies[type_of_module_pair[1]])
        eff *= (detection_bin_efficiencies0[detection_bin_1] *
                detection_bin_efficiencies1[detection_bin_2])
        if eff == 0:
            return 0.

    # per module-pair efficiencies
    module_pair_efficiencies_vectors = (
        scanner.detection_efficiencies.module_pair_efficiencies_vectors)
    if module_pair_efficiencies_vectors is not None:
        module_pair_SGID_LUT = scanner.detection_efficiencies.module_pair_sgidlut
        assert module_pair_SGID_LUT is not None
        expanded_det_bin0 = expand_detection_bin(scanner,
                                                 type_of_module_pair[0],
                                                 detection_bin_1)
        expanded_det_bin1 = expand_detection_bin(scanner,
                                                 type_of_module_pair[1],
                                                 detection_bin_2)

        SGID = module_pair_SGID_LUT[type_of_module_pair[0]][
            type_of_module_pair[1]][expanded_det_bin0.module_index,
                                    expanded_det_bin1.module_index]
        if SGID < 0:
            return 0.
        module_pair_efficiencies = module_pair_efficiencies_vectors[
            type_of_module_pair[0]][type_of_module_pair[1]][SGID]
        assert module_pair_efficiencies.sgid == SGID
        num_en0 = scanner.event_energy_bin_edges[
            type_of_module_pair[0]].number_of_bins()
        num_en1 = scanner.event_energy_bin_edges[
            type_of_module_pair[1]].number_of_bins()
        # TODO create helper for next calculation
        eff *= module_pair_efficiencies.values[
            expanded_det_bin0.element_index * num_en0 +
            expanded_det_bin0.energy_index,
            expanded_det_bin1.element_index * num_en1 +
            expanded_det_bin1.energy_index]

    return eff
