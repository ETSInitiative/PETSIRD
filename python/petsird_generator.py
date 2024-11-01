#  Copyright (C) 2022-2023 Microsoft Corporation
#  Copyright (C) 2023-2024 University College London
#
#  SPDX-License-Identifier: Apache-2.0

import sys
import math
import random
import numpy
import numpy.typing as npt
from collections.abc import Iterator
import petsird
from petsird_helpers import get_num_det_els, get_module_and_element

# these are constants for now
NUMBER_OF_ENERGY_BINS = 3
NUMBER_OF_TOF_BINS = 300
RADIUS = 400
CRYSTAL_LENGTH = (20, 4, 4)
# num crystals in a module
NUM_CRYSTALS_PER_MODULE = (2, 4,5)
NUM_MODULES_ALONG_RING = 20
NUM_MODULES_ALONG_AXIS = 2
MODULE_AXIS_SPACING = (NUM_CRYSTALS_PER_MODULE[2] + 4) * CRYSTAL_LENGTH[2]
NUMBER_OF_TIME_BLOCKS = 6
NUMBER_OF_EVENTS = 1000
COUNT_RATE = 500


def make_coordinate(v: tuple) -> petsird.Coordinate:
    return petsird.Coordinate(c=numpy.array(v, dtype=numpy.float32))


def get_crystal() -> petsird.SolidVolume:
    """return a cuboid volume"""
    crystal_shape = petsird.BoxShape(
        corners=[
            make_coordinate((0, 0, 0)),
            make_coordinate((0, 0, CRYSTAL_LENGTH[2])),
            make_coordinate((0, CRYSTAL_LENGTH[1], CRYSTAL_LENGTH[2])),
            make_coordinate((0, CRYSTAL_LENGTH[1], 0)),
            make_coordinate((CRYSTAL_LENGTH[0], 0, 0)),
            make_coordinate((CRYSTAL_LENGTH[0], 0, CRYSTAL_LENGTH[2])),
            make_coordinate((CRYSTAL_LENGTH[0], CRYSTAL_LENGTH[1], CRYSTAL_LENGTH[2])),
            make_coordinate((CRYSTAL_LENGTH[0], CRYSTAL_LENGTH[1], 0)),
        ]
    )
    return petsird.BoxSolidVolume(shape=crystal_shape, material_id=1)


def get_detector_module() -> petsird.DetectorModule:
    """return a module of NUM_CRYSTALS_PER_MODULE cuboids"""
    crystal = get_crystal()
    rep_volume = petsird.ReplicatedBoxSolidVolume(object=crystal)
    N0 = NUM_CRYSTALS_PER_MODULE[0]
    N1 = NUM_CRYSTALS_PER_MODULE[1]
    N2 = NUM_CRYSTALS_PER_MODULE[2]
    for rep0 in range(N0):
        for rep1 in range(N1):
            for rep2 in range(N2):
                transform = petsird.RigidTransformation(
                    matrix=numpy.array(
                        (
                            (1.0, 0.0, 0.0, RADIUS + rep0 * CRYSTAL_LENGTH[0]),
                            (0.0, 1.0, 0.0, (rep1 - N1/2) * CRYSTAL_LENGTH[1]),
                            (0.0, 0.0, 1.0, (rep2 - N2/2) * CRYSTAL_LENGTH[2]),
                        ),
                        dtype="float32",
                    )
                )
                rep_volume.transforms.append(transform)
                rep_volume.ids.append(rep0 + N0 * (rep1 + N1 * rep2))

    return petsird.DetectorModule(
        detecting_elements=[rep_volume], detecting_element_ids=[0]
    )


def get_scanner_geometry() -> petsird.ScannerGeometry:
    """return a scanner build by rotating a module around the (0,0,1) axis"""
    detector_module = get_detector_module()
    radius = RADIUS
    angles = [2 * math.pi * i / NUM_MODULES_ALONG_RING for i in range(NUM_MODULES_ALONG_RING)]

    rep_module = petsird.ReplicatedDetectorModule(object=detector_module)
    module_id = 0
    for angle in angles:
        for ax_mod in range(NUM_MODULES_ALONG_AXIS):
            transform = petsird.RigidTransformation(
                matrix=numpy.array(
                    (
                        (math.cos(angle), math.sin(angle), 0.0, 0.0),
                        (-math.sin(angle), math.cos(angle), 0.0, 0.0),
                        (0.0, 0.0, 1.0, MODULE_AXIS_SPACING*ax_mod),
                    ),
                    dtype="float32",
                )
            )
            rep_module.ids.append(module_id)
            module_id += 1
            rep_module.transforms.append(transform)

    return petsird.ScannerGeometry(replicated_modules=[rep_module], ids=[0])


def get_detection_efficiencies(
    scanner: petsird.ScannerInformation,
) -> petsird.DetectionEfficiencies:
    num_det_els = get_num_det_els(scanner.scanner_geometry)
    det_el_efficiencies = numpy.ones(
        (num_det_els, scanner.number_of_energy_bins()), dtype=numpy.float32
    )

    # only 1 type of module in the current scanner
    assert len(scanner.scanner_geometry.replicated_modules) == 1
    rep_module = scanner.scanner_geometry.replicated_modules[0]
    num_modules = len(rep_module.transforms)
    # We will only use rotational symmetries (no translation along the axis yet)
    # We also assume all module-pairs are in coincidence, except those with the same angle.
    # Writing a module number as (z-position, angle):
    #   eff((z1,a1), (z2, a2)) == eff((z1,0), (z2, abs(a2-a1)))
    # or in linear indices
    #   eff(z1 + NZ * a1, z2 + NZ * a2) == eff(z1, z2 + NZ * abs(a2 - a1))
    # (coincident) SGIDs need to start from 0, so ignoring self-coincident angles
    num_SGIDs = NUM_MODULES_ALONG_AXIS * NUM_MODULES_ALONG_AXIS * (NUM_MODULES_ALONG_RING-1)
    # SGID = z1 + NZ * (z2 + NZ * abs(a2 - a1) - 1)
    NZ = NUM_MODULES_ALONG_AXIS
    module_pair_SGID_LUT = numpy.ndarray((num_modules, num_modules), dtype="int32")
    for mod1 in range(num_modules):
        for mod2 in range(num_modules):
            z1 = mod1 % NZ
            a1 = mod1 // NZ
            z2 = mod2 % NZ
            a2 = mod2 // NZ
            if a1 == a2:
                module_pair_SGID_LUT[mod1, mod2]= -1
            else:
                module_pair_SGID_LUT[mod1, mod2] = z1 + NZ * (z2 + NZ * (abs(a2 - a1) - 1))

    # print("SGID LUT:\n", module_pair_SGID_LUT, file=sys.stderr)
    assert numpy.max(module_pair_SGID_LUT) == num_SGIDs - 1
    module_pair_efficiencies_vector = []
    assert len(rep_module.object.detecting_elements) == 1
    detecting_elements = rep_module.object.detecting_elements[0]
    num_det_els_in_module = len(detecting_elements.transforms)
    for SGID in range(num_SGIDs):
        # extract first module_pair for this SGID. However, as this currently unused, it is commented out
        # module_pair = numpy.argwhere(module_pair_SGID_LUT == SGID)[0]
        # print(module_pair, file=sys.stderr)
        module_pair_efficiencies = numpy.ones(
            (
                num_det_els_in_module,
                scanner.number_of_energy_bins(),
                num_det_els_in_module,
                scanner.number_of_energy_bins(),
            ),
            dtype=numpy.float32,
        )
        # give some (non-physical) value
        module_pair_efficiencies *= SGID
        module_pair_efficiencies_vector.append(
            petsird.ModulePairEfficiencies(values=module_pair_efficiencies, sgid=SGID)
        )
        assert len(module_pair_efficiencies_vector) == SGID + 1

    return petsird.DetectionEfficiencies(
        det_el_efficiencies=det_el_efficiencies,
        module_pair_sgidlut=module_pair_SGID_LUT,
        module_pair_efficiencies_vector=module_pair_efficiencies_vector,
    )


def get_scanner_info() -> petsird.ScannerInformation:

    scanner_geometry = get_scanner_geometry()

    # TODO scanner_info.bulk_materials

    # TOF info (in mm)
    tofBinEdges = numpy.linspace(
        -RADIUS, RADIUS, NUMBER_OF_TOF_BINS + 1, dtype="float32"
    )
    energyBinEdges = numpy.linspace(
        430, 650, NUMBER_OF_ENERGY_BINS + 1, dtype="float32"
    )
    # need energy bin info before being able to construct the detection efficiencies
    # so we will construct a scanner without the efficiencies first
    scanner = petsird.ScannerInformation(
        model_name="PETSIRD_TEST",
        scanner_geometry=scanner_geometry,
        tof_bin_edges=tofBinEdges,
        tof_resolution=9.4,  # in mm
        energy_bin_edges=energyBinEdges,
        energy_resolution_at_511=0.11,  # as fraction of 511
        event_time_block_duration=1,  # ms
    )

    scanner.detection_efficiencies = get_detection_efficiencies(scanner)
    return scanner


def get_header() -> petsird.Header:
    subject = petsird.Subject(id="123456")
    institution = petsird.Institution(
        name="Diamond Light Source",
        address="Harwell Science and Innovation Campus, Didcot, Oxfordshire, OX11 0DE, UK",
    )
    return petsird.Header(
        exam=petsird.ExamInformation(subject=subject, institution=institution),
        scanner=get_scanner_info(),
    )


def get_events(
    header: petsird.Header, num_events: int
) -> Iterator[petsird.CoincidenceEvent]:
    """Generate some random events"""
    detector_count = get_num_det_els(header.scanner.scanner_geometry)
    for _ in range(num_events):
        # Generate random detector_ids, where the corresponding modules are in coincidence
        while True:
            detector_ids = [
                random.randrange(0, detector_count),
                random.randrange(0, detector_count),
            ]
            mod_and_els = get_module_and_element(
                header.scanner.scanner_geometry, detector_ids
            )
            if header.scanner.detection_efficiencies.module_pair_sgidlut[mod_and_els[0].module, mod_and_els[1].module] >= 0:
                # in coincidence, we can get out of the loop
                break

        yield petsird.CoincidenceEvent(
            detector_ids=detector_ids,
            energy_indices=[
                random.randrange(0, NUMBER_OF_ENERGY_BINS),
                random.randrange(0, NUMBER_OF_ENERGY_BINS),
            ],
            tof_idx=random.randrange(0, NUMBER_OF_TOF_BINS),
        )


if __name__ == "__main__":
    # numpy random number generator
    rng = numpy.random.default_rng()

    with petsird.BinaryPETSIRDWriter(sys.stdout.buffer) as writer:
        # with petsird.NDJsonPETSIRDWriter(sys.stdout) as writer:
        header = get_header()
        writer.write_header(header)
        for t in range(NUMBER_OF_TIME_BLOCKS):
            start = t * header.scanner.event_time_block_duration
            num_prompts_this_block = rng.poisson(COUNT_RATE)
            prompts_this_block = list(get_events(header, num_prompts_this_block))
            # Normally we'd write multiple blocks, but here we have just one, so let's write a tuple with just one element
            writer.write_time_blocks(
                (
                    petsird.TimeBlock.EventTimeBlock(
                        petsird.EventTimeBlock(
                            start=start, prompt_events=prompts_this_block
                        )
                    ),
                )
            )
