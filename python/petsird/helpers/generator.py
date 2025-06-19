#  Copyright (C) 2022-2023 Microsoft Corporation
#  Copyright (C) 2023-2025 University College London
#
#  SPDX-License-Identifier: Apache-2.0
import math
import random
import sys
from collections.abc import Iterator

import numpy
import petsird
from petsird.helpers import get_detection_efficiency, get_num_detection_bins

# these are constants for now
NUMBER_OF_EVENT_ENERGY_BINS = 3
NUMBER_OF_TOF_BINS = 300
RADIUS = 400
CRYSTAL_LENGTH = (20, 4, 4)
# num crystals in a module
NUM_CRYSTALS_PER_MODULE = (2, 4, 7)
NUM_MODULES_ALONG_RING = 20
NUM_MODULES_ALONG_AXIS = 2
MODULE_AXIS_SPACING = (NUM_CRYSTALS_PER_MODULE[2] + 4) * CRYSTAL_LENGTH[2]
NUMBER_OF_TIME_BLOCKS = 6
NUMBER_OF_EVENTS = 1000
COUNT_RATE = 500  # 1/ms
EVENT_TIME_BLOCK_DURATION = 1  # ms


def make_coordinate(v: tuple) -> petsird.Coordinate:
    return petsird.Coordinate(c=numpy.array(v, dtype=numpy.float32))


def get_crystal() -> petsird.SolidVolume:
    """return a cuboid volume"""
    crystal_shape = petsird.BoxShape(corners=[
        make_coordinate((0, 0, 0)),
        make_coordinate((0, 0, CRYSTAL_LENGTH[2])),
        make_coordinate((0, CRYSTAL_LENGTH[1], CRYSTAL_LENGTH[2])),
        make_coordinate((0, CRYSTAL_LENGTH[1], 0)),
        make_coordinate((CRYSTAL_LENGTH[0], 0, 0)),
        make_coordinate((CRYSTAL_LENGTH[0], 0, CRYSTAL_LENGTH[2])),
        make_coordinate((CRYSTAL_LENGTH[0], CRYSTAL_LENGTH[1],
                         CRYSTAL_LENGTH[2])),
        make_coordinate((CRYSTAL_LENGTH[0], CRYSTAL_LENGTH[1], 0)),
    ])
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
                transform = petsird.RigidTransformation(matrix=numpy.array(
                    (
                        (1.0, 0.0, 0.0, RADIUS + rep0 * CRYSTAL_LENGTH[0]),
                        (0.0, 1.0, 0.0, (rep1 - N1 / 2) * CRYSTAL_LENGTH[1]),
                        (0.0, 0.0, 1.0, (rep2 - N2 / 2) * CRYSTAL_LENGTH[2]),
                    ),
                    dtype="float32",
                ))
                rep_volume.transforms.append(transform)

    return petsird.DetectorModule(detecting_elements=rep_volume)


def get_scanner_geometry() -> petsird.ScannerGeometry:
    """return a scanner build by rotating a module around the (0,0,1) axis"""
    detector_module = get_detector_module()
    angles = [
        2 * math.pi * i / NUM_MODULES_ALONG_RING
        for i in range(NUM_MODULES_ALONG_RING)
    ]

    rep_module = petsird.ReplicatedDetectorModule(object=detector_module)
    for angle in angles:
        for ax_mod in range(NUM_MODULES_ALONG_AXIS):
            transform = petsird.RigidTransformation(matrix=numpy.array(
                (
                    (math.cos(angle), math.sin(angle), 0.0, 0.0),
                    (-math.sin(angle), math.cos(angle), 0.0, 0.0),
                    (0.0, 0.0, 1.0, MODULE_AXIS_SPACING * ax_mod),
                ),
                dtype="float32",
            ))
            rep_module.transforms.append(transform)

    return petsird.ScannerGeometry(replicated_modules=[rep_module])


def get_detection_efficiencies(
    scanner: petsird.ScannerInformation, ) -> petsird.DetectionEfficiencies:
    """return some (non-physical) detection efficiencies"""
    # TODO only 1 type of module in the current scanner
    assert len(scanner.scanner_geometry.replicated_modules) == 1
    type_of_module = 0
    num_detection_bins = get_num_detection_bins(scanner, type_of_module)
    event_energy_bin_edges = scanner.event_energy_bin_edges[type_of_module]
    num_event_energy_bins = event_energy_bin_edges.number_of_bins()
    detection_bin_efficiencies = numpy.ones((num_detection_bins),
                                            dtype=numpy.float32)

    rep_module = scanner.scanner_geometry.replicated_modules[type_of_module]
    num_modules = len(rep_module.transforms)
    # We will use rotational symmetries translation along the axis
    # We assume all module-pairs are in coincidence, except those
    # with the same angle.
    # Writing a module number as (z-position, angle):
    #   eff((z1,a1), (z2, a2)) == eff((z1,0), (z2, abs(a2-a1)))
    # or in linear indices
    #   eff(z1 + NZ * a1, z2 + NZ * a2) == eff(z1, z2 + NZ * abs(a2 - a1))
    # (coincident) SGIDs need to start from 0, so ignoring self-coincident
    # angles we get
    #   SGID = z1 + NZ * (z2 + NZ * abs(a2 - a1) - 1)
    num_SGIDs = NUM_MODULES_ALONG_AXIS * NUM_MODULES_ALONG_AXIS * (
        NUM_MODULES_ALONG_RING - 1)
    NZ = NUM_MODULES_ALONG_AXIS
    module_pair_SGID_LUT = numpy.ndarray((num_modules, num_modules),
                                         dtype="int32")
    for mod1 in range(num_modules):
        for mod2 in range(num_modules):
            z1 = mod1 % NZ
            a1 = mod1 // NZ
            z2 = mod2 % NZ
            a2 = mod2 // NZ
            if a1 == a2:
                module_pair_SGID_LUT[mod1, mod2] = -1
            else:
                module_pair_SGID_LUT[mod1,
                                     mod2] = z1 + NZ * (z2 + NZ *
                                                        (abs(a2 - a1) - 1))

    # print("SGID LUT:\n", module_pair_SGID_LUT, file=sys.stderr)
    assert numpy.max(module_pair_SGID_LUT) == num_SGIDs - 1
    module_pair_efficiencies_vector = []
    detecting_elements = rep_module.object.detecting_elements
    num_det_els_in_module = len(detecting_elements.transforms)
    event_energy_bin_edges = scanner.event_energy_bin_edges[type_of_module]
    num_event_energy_bins = event_energy_bin_edges.number_of_bins()

    for SGID in range(num_SGIDs):
        # Extract first module_pair for this SGID. However, as this
        # currently unused, it is commented out
        # module_pair = numpy.argwhere(module_pair_SGID_LUT == SGID)[0]
        # print(module_pair, file=sys.stderr)
        module_pair_efficiencies = numpy.ones(
            (
                num_det_els_in_module * num_event_energy_bins,
                num_det_els_in_module * num_event_energy_bins,
            ),
            dtype=numpy.float32,
        )
        # give some (non-physical) value
        module_pair_efficiencies *= SGID
        module_pair_efficiencies_vector.append(
            petsird.ModulePairEfficiencies(values=module_pair_efficiencies,
                                           sgid=SGID))
        assert len(module_pair_efficiencies_vector) == SGID + 1

    return petsird.DetectionEfficiencies(
        detection_bin_efficiencies=[detection_bin_efficiencies],
        module_pair_sgidlut=[[module_pair_SGID_LUT]],
        module_pair_efficiencies_vectors=[[module_pair_efficiencies_vector]],
    )


def get_scanner_info() -> petsird.ScannerInformation:

    scanner_geometry = get_scanner_geometry()
    num_types_of_modules = scanner_geometry.number_of_replicated_modules()
    # Example code below is restricted to the following case
    assert num_types_of_modules == 1

    # TODO scanner_info.bulk_materials

    # TOF info (in mm)
    tofBinEdges = petsird.BinEdges(edges=numpy.linspace(
        -RADIUS, RADIUS, NUMBER_OF_TOF_BINS + 1, dtype="float32"))
    energyBinEdges = petsird.BinEdges(edges=numpy.linspace(
        430, 650, NUMBER_OF_EVENT_ENERGY_BINS + 1, dtype="float32"))
    # In this example, there is only 1 module-type, therefore we need 1x1 "nested lists"
    allTofBinEdges = [[tofBinEdges]]
    tofResolution = [[9.4]]  # in mm
    # 1-element list for energy info
    allEnergyBinEdges = [energyBinEdges]
    energyResolutionAt511 = [0.11]  # as fraction of 511

    # We need energy bin info before being able to construct the detection
    # efficiencies, so we first construct a scanner without the efficiencies
    scanner = petsird.ScannerInformation(
        model_name="PETSIRD_TEST",
        scanner_geometry=scanner_geometry,
        tof_bin_edges=allTofBinEdges,
        tof_resolution=tofResolution,
        event_energy_bin_edges=allEnergyBinEdges,
        energy_resolution_at_511=energyResolutionAt511)

    # Now added the efficiencies
    scanner.detection_efficiencies = get_detection_efficiencies(scanner)

    scanner.coincidence_policy = petsird.CoincidencePolicy.REJECT_MULTIPLES
    scanner.delayed_coincidences_are_stored = False
    scanner.triple_events_are_stored = False

    return scanner


def get_header() -> petsird.Header:
    subject = petsird.Subject(id="123456")
    institution = petsird.Institution(
        name="Some institution",
        address="Didcot, Oxfordshire, OX11 0DE, UK",
    )
    return petsird.Header(
        exam=petsird.ExamInformation(subject=subject, institution=institution),
        scanner=get_scanner_info(),
    )


def get_random_uint(max):
    return random.randrange(0, max)


def get_events(header: petsird.Header,
               type_of_module_pair: petsird.TypeOfModulePair,
               num_events: int) -> Iterator[petsird.CoincidenceEvent]:
    """Generate some random events"""
    num_modules0 = len(header.scanner.scanner_geometry.replicated_modules[
        type_of_module_pair[0]].transforms)
    num_detecting_elements0 = len(
        header.scanner.scanner_geometry.replicated_modules[
            type_of_module_pair[0]].object.detecting_elements.transforms)
    event_energy_bin_edges0 = header.scanner.event_energy_bin_edges[
        type_of_module_pair[0]]
    num_event_energy_bins0 = event_energy_bin_edges0.number_of_bins()
    tof_bin_edges = header.scanner.tof_bin_edges[type_of_module_pair[0]][
        type_of_module_pair[1]]
    num_tof_bins = tof_bin_edges.number_of_bins()
    count1 = get_num_detection_bins(header.scanner, type_of_module_pair[1])

    detection_bins = [petsird.DetectionBin(), petsird.DetectionBin()]
    for _ in range(num_events):
        event = petsird.CoincidenceEvent(
            detection_bins=detection_bins,
            tof_idx=get_random_uint(num_tof_bins),
        )
        # Generate random detection_bins until detection effficiency is not zero
        while True:
            expanded_detection_bin0 = petsird.ExpandedDetectionBin(
                module_index=get_random_uint(num_modules0),
                element_index=get_random_uint(num_detecting_elements0),
                energy_index=get_random_uint(num_event_energy_bins0),
            )
            event.detection_bins[0] = petsird.helpers.make_detection_bin(
                header.scanner, type_of_module_pair[0],
                expanded_detection_bin0)
            # TODO move test to separate function
            assert expanded_detection_bin0 == petsird.helpers.expand_detection_bin(
                header.scanner, type_of_module_pair[0], detection_bins[0])

            # short-cut to directly generate a random detection bin
            event.detection_bins[1] = get_random_uint(count1)

            if get_detection_efficiency(header.scanner, type_of_module_pair,
                                        event) > 0:
                # in coincidence, we can get out of the loop
                break

        yield event


if __name__ == "__main__":
    # numpy random number generator
    rng = numpy.random.default_rng()

    with petsird.BinaryPETSIRDWriter(sys.stdout.buffer) as writer:
        # with petsird.NDJsonPETSIRDWriter(sys.stdout) as writer:
        header = get_header()
        writer.write_header(header)
        for t in range(NUMBER_OF_TIME_BLOCKS):
            time_interval = petsird.TimeInterval(
                start=t * EVENT_TIME_BLOCK_DURATION,
                stop=(t + 1) * EVENT_TIME_BLOCK_DURATION)
            average_num = EVENT_TIME_BLOCK_DURATION * COUNT_RATE
            num_prompts_this_block = rng.poisson(average_num)
            type_of_module_pair = petsird.TypeOfModulePair((0, 0))
            prompts_this_block = [[
                list(
                    get_events(header, type_of_module_pair,
                               num_prompts_this_block))
            ]]
            # Normally we'd write multiple blocks, but here we have just one,
            # so let's write a tuple with just one element
            writer.write_time_blocks((petsird.TimeBlock.EventTimeBlock(
                petsird.EventTimeBlock(time_interval=time_interval,
                                       prompt_events=prompts_this_block)), ))
