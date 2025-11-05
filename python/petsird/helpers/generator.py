#  Copyright (C) 2022-2023 Microsoft Corporation
#  Copyright (C) 2023-2025 University College London
#
#  SPDX-License-Identifier: Apache-2.0
"""This is an example file how to create a PETSIRD file.

It creates a demo scanner with 2 types of modules, and some random events.
Efficiencies are also stores, although the values are meaningless.

This code only serves as illustration on how to create a PETSIRD file, and
will need serious adaption to be useful.
"""
import math
import random
import sys
from collections.abc import Iterator
from dataclasses import dataclass
from typing import Iterable, List, Tuple

import numpy
import petsird
from petsird.helpers import get_detection_efficiency, get_num_detection_bins


@dataclass
class CylindricalBlocksInfo:
    """ Definitions for a cylindrical arrangement of modules of one particular type.

    This is incomplete of course, but hopefully sufficient for illustration.

    Definitions here also include TOF, which is supposed for "this module-type"
    coincidences.
    """
    crystal_length: Tuple[float]
    num_crystals_per_module: Tuple[int]
    num_modules_along_ring: int
    num_modules_along_axis: int
    radius: float  # in mm
    arc: float  # in radians
    start_angle: float  # in radians
    module_spacing_along_axis: float
    number_of_tof_bins: int
    tof_resolution: float  # CTR in mm for coincidences between modules of this type
    LLD: float  # lower level discriminator (keV)
    ULD: float  # upper level discriminator (keV)
    number_of_event_energy_bins: int
    energy_resolution: float  # FWHM at 511 keV as fraction of 511.
    material_id: int


# 2 examples, not very physical plausible/efficient.
# We will build a scanner from both below.
mtype0_def = CylindricalBlocksInfo(
    radius=400,
    crystal_length=(20, 4, 4),
    num_crystals_per_module=(2, 4, 7),
    num_modules_along_ring=20,
    num_modules_along_axis=2,
    arc=2 * math.pi,
    start_angle=0,
    module_spacing_along_axis=(7 + 24) *
    4.,  # (num_crystals_per_module[2] + 24) * crystal_length[2],
    number_of_tof_bins=11,
    tof_resolution=9,
    LLD=450,
    ULD=650,
    number_of_event_energy_bins=3,
    energy_resolution=0.1,
    material_id=1,
)
mtype1_def = CylindricalBlocksInfo(
    radius=150,
    crystal_length=(10, 2, 2),
    num_crystals_per_module=(1, 3, 30),
    num_modules_along_ring=15,
    num_modules_along_axis=1,
    arc=math.pi,
    start_angle=0,
    module_spacing_along_axis=0.,
    number_of_tof_bins=5,
    tof_resolution=6,
    LLD=460,
    ULD=640,
    number_of_event_energy_bins=1,
    energy_resolution=0.8,
    material_id=1,
)
# Some constants related to how many events we will generate etc
NUMBER_OF_TIME_BLOCKS = 6
NUMBER_OF_EVENTS = 1000
COUNT_RATE = 500  # 1/ms
EVENT_TIME_BLOCK_DURATION = 1  # ms


def make_coordinate(v: tuple) -> petsird.Coordinate:
    return petsird.Coordinate(c=numpy.array(v, dtype=numpy.float32))


def get_crystal(crystal_length: Tuple[float],
                material_id: int = 1) -> petsird.SolidVolume:
    """return a cuboid volume with first corner at 0,0,0"""
    crystal_shape = petsird.BoxShape(corners=[
        make_coordinate((0, 0, 0)),
        make_coordinate((0, 0, crystal_length[2])),
        make_coordinate((0, crystal_length[1], crystal_length[2])),
        make_coordinate((0, crystal_length[1], 0)),
        make_coordinate((crystal_length[0], 0, 0)),
        make_coordinate((crystal_length[0], 0, crystal_length[2])),
        make_coordinate((crystal_length[0], crystal_length[1],
                         crystal_length[2])),
        make_coordinate((crystal_length[0], crystal_length[1], 0)),
    ])
    return petsird.BoxSolidVolume(shape=crystal_shape, material_id=material_id)


def get_detector_module(
        module_def: CylindricalBlocksInfo) -> petsird.DetectorModule:
    r"""return a module of num_crystals_per_module cuboids

    Currently transform are such that the first axis is the "crystal length",
    radius is the "inner radius", and the crystal is centered along
    the other directions.
    """
    crystal = get_crystal(module_def.crystal_length, module_def.material_id)
    rep_volume = petsird.ReplicatedBoxSolidVolume(object=crystal)
    N0 = module_def.num_crystals_per_module[0]
    N1 = module_def.num_crystals_per_module[1]
    N2 = module_def.num_crystals_per_module[2]
    for rep0 in range(N0):
        for rep1 in range(N1):
            for rep2 in range(N2):
                transform = petsird.RigidTransformation(matrix=numpy.array(
                    (
                        (1.0, 0.0, 0.0, module_def.radius +
                         rep0 * module_def.crystal_length[0]),
                        (0.0, 1.0, 0.0,
                         (rep1 - (N1 - 1) / 2) * module_def.crystal_length[1]),
                        (0.0, 0.0, 1.0,
                         (rep2 - (N2 - 1) / 2) * module_def.crystal_length[2]),
                    ),
                    dtype="float32",
                ))
                rep_volume.transforms.append(transform)

    return petsird.DetectorModule(detecting_elements=rep_volume)


def get_replicated_module(
        module_def: CylindricalBlocksInfo) -> petsird.ReplicatedDetectorModule:
    """rotating a module around, translating along the (0,0,1) axis

    The "centre" of the translated/rotated modules is at (0,0,0) (if arc=2pi)."""
    detector_module = get_detector_module(module_def)
    angles = [
        module_def.start_angle +
        module_def.arc * i / module_def.num_modules_along_ring
        for i in range(module_def.num_modules_along_ring)
    ]

    N2 = module_def.num_modules_along_axis
    rep_module = petsird.ReplicatedDetectorModule(object=detector_module)
    for angle in angles:
        for ax_mod in range(N2):
            transform = petsird.RigidTransformation(matrix=numpy.array(
                (
                    (math.cos(angle), math.sin(angle), 0.0, 0.0),
                    (-math.sin(angle), math.cos(angle), 0.0, 0.0),
                    (0.0, 0.0, 1.0, module_def.module_spacing_along_axis *
                     (ax_mod - (N2 - 1) / 2)),
                ),
                dtype="float32",
            ))
            rep_module.transforms.append(transform)
    return rep_module


def get_scanner_geometry(
        module_defs: Iterable[CylindricalBlocksInfo]
) -> petsird.ScannerGeometry:
    """return a scanner build by rotating one or more modules around the (0,0,1) axis"""

    return petsird.ScannerGeometry(replicated_modules=[
        get_replicated_module(m_def) for m_def in module_defs
    ])


def get_detection_bin_efficiencies(
    scanner: petsird.ScannerInformation,
    type_of_module: int,
    module_def: CylindricalBlocksInfo,
) -> numpy.ndarray:
    """return some (non-physical) detection bin efficiencies"""
    num_detection_bins = get_num_detection_bins(scanner, type_of_module)
    detection_bin_efficiencies = numpy.ones((num_detection_bins),
                                            dtype=numpy.float32)
    return detection_bin_efficiencies


def create_empty_module_pair_SGID_LUT(scanner: petsird.ScannerInformation,
                                      type_of_module0: int,
                                      type_of_module1: int) -> numpy.ndarray:
    """return a LUT that says no modules are in coincidence"""
    num_modules0 = len(scanner.scanner_geometry.
                       replicated_modules[type_of_module0].transforms)
    num_modules1 = len(scanner.scanner_geometry.
                       replicated_modules[type_of_module1].transforms)
    LUT = numpy.ndarray((num_modules0, num_modules1), dtype="int32")
    LUT[:] = -1
    return LUT


def create_empty_module_pair_efficiencies(
        scanner: petsird.ScannerInformation, type_of_module0: int,
        type_of_module1: int) -> numpy.ndarray:
    """return a 2d array (filled with 1) of the appropriate size"""
    rep_module = scanner.scanner_geometry.replicated_modules[type_of_module0]
    detecting_elements0 = rep_module.object.detecting_elements
    num_det_els_in_module0 = len(detecting_elements0.transforms)
    event_energy_bin_edges0 = scanner.event_energy_bin_edges[type_of_module0]
    num_event_energy_bins0 = event_energy_bin_edges0.number_of_bins()
    rep_module = scanner.scanner_geometry.replicated_modules[type_of_module1]
    detecting_elements1 = rep_module.object.detecting_elements
    num_det_els_in_module1 = len(detecting_elements1.transforms)
    event_energy_bin_edges1 = scanner.event_energy_bin_edges[type_of_module1]
    num_event_energy_bins1 = event_energy_bin_edges1.number_of_bins()
    module_pair_efficiencies = numpy.ones(
        (
            num_det_els_in_module0 * num_event_energy_bins0,
            num_det_els_in_module1 * num_event_energy_bins1,
        ),
        dtype=numpy.float32,
    )
    return module_pair_efficiencies


def get_module_pair_efficiencies_one_module_type(
    scanner: petsird.ScannerInformation,
    type_of_module: int,
    module_def: CylindricalBlocksInfo,
) -> Tuple[numpy.ndarray, List[petsird.ModulePairEfficiencies]]:
    """return detection efficiencies for a module-pair of the same type

    The function returns a tuple with module_pair_SGID_LUT,
    module_pair_efficiencies_vector.

    At present, all detector modules are assumed in coincidence (except with itself).
    Actual values are non-sensical.
    """
    rep_module = scanner.scanner_geometry.replicated_modules[type_of_module]
    num_modules = len(rep_module.transforms)

    # Ideally we would get the following from the transforms, but we need to
    # know the scanner geometry anyway for the symmetries, so we get it from module_def
    num_modules_along_axis = module_def.num_modules_along_axis
    num_modules_along_ring = module_def.num_modules_along_ring
    # We will (only) use rotational symmetries translation along the axis
    # We assume all module-pairs are in coincidence, except those
    # with the same angle.
    # Writing a module number as (z-position, angle):
    #   eff((z1,a1), (z2, a2)) == eff((z1,0), (z2, abs(a2-a1)))
    # or in linear indices
    #   eff(z1 + NZ * a1, z2 + NZ * a2) == eff(z1, z2 + NZ * abs(a2 - a1))
    # (coincident) SGIDs need to start from 0, so ignoring self-coincident
    # angles we get
    #   SGID = z1 + NZ * (z2 + NZ * abs(a2 - a1) - 1)
    num_SGIDs = num_modules_along_axis * num_modules_along_axis * (
        num_modules_along_ring - 1)
    NZ = num_modules_along_axis
    module_pair_SGID_LUT = create_empty_module_pair_SGID_LUT(
        scanner, type_of_module, type_of_module)
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

    module_pair_efficiencies = create_empty_module_pair_efficiencies(
        scanner, type_of_module, type_of_module)
    for SGID in range(num_SGIDs):
        # Extract first module_pair for this SGID. However, as this is
        # currently unused, it is commented out
        # module_pair = numpy.argwhere(module_pair_SGID_LUT == SGID)[0]
        # print(module_pair, file=sys.stderr)
        # give some (non-physical) value
        module_pair_efficiencies_vector.append(
            petsird.ModulePairEfficiencies(values=module_pair_efficiencies *
                                           SGID,
                                           sgid=SGID))
        assert len(module_pair_efficiencies_vector) == SGID + 1

    return (module_pair_SGID_LUT, module_pair_efficiencies_vector)


def get_module_pair_efficiencies_two_module_types(
    scanner: petsird.ScannerInformation,
    type_of_module0: int,
    type_of_module1: int,
    module_def0: CylindricalBlocksInfo,
    module_def1: CylindricalBlocksInfo,
) -> Tuple[numpy.ndarray, List[petsird.ModulePairEfficiencies]]:
    """return detection efficiencies for a module-pair of different types

    The function returns a tuple with module_pair_SGID_LUT,
    module_pair_efficiencies_vector.

    As discovering symmetries here is hard, we'll just ignore them.

    Actual values are non-sensical.
    """
    LUT = create_empty_module_pair_SGID_LUT(scanner, type_of_module0,
                                            type_of_module1)
    module_pair_efficiencies = create_empty_module_pair_efficiencies(
        scanner, type_of_module0, type_of_module1)
    SGID = 0  # first SGID needs to be 0
    module_pair_efficiencies_vector = []
    for m0 in range(LUT.shape[0]):
        for m1 in range(LUT.shape[1]):
            LUT[m0, m1] = SGID
            # using nonsensical values for the following
            petsird_module_pair_efficiencies = petsird.ModulePairEfficiencies(
                values=module_pair_efficiencies * SGID, sgid=SGID)
            module_pair_efficiencies_vector.append(
                petsird_module_pair_efficiencies)
            SGID += 1

    return LUT, module_pair_efficiencies_vector


def get_detection_efficiencies(
    scanner: petsird.ScannerInformation,
    module_defs: Iterable[CylindricalBlocksInfo],
) -> petsird.DetectionEfficiencies:
    """return some (non-physical) detection efficiencies"""

    # use some non-physical value for the calibration factor
    calibration_factor = 42.
    num_types_of_modules = scanner.scanner_geometry.number_of_module_types()
    # single-level list for detection_bin_efficiencies
    all_detection_bin_efficiencies = []
    # double-level list for module_pair things
    all_LUTs = []
    all_module_pair_efficiencies_vectors = []
    for mtype0 in range(num_types_of_modules):
        detection_bin_efficiencies = get_detection_bin_efficiencies(
            scanner, mtype0, module_defs[mtype0])
        all_detection_bin_efficiencies.append(detection_bin_efficiencies)
        # single-level lists for the information at the current mtype0
        LUTs0 = []
        module_pair_efficiencies_vectors0 = []
        for mtype1 in range(num_types_of_modules):
            # construct information
            if (mtype0 == mtype1):
                (LUT, module_pair_efficiencies_vector
                 ) = get_module_pair_efficiencies_one_module_type(
                     scanner, mtype0, module_defs[mtype0])
            else:
                (LUT, module_pair_efficiencies_vector
                 ) = get_module_pair_efficiencies_two_module_types(
                     scanner, mtype0, mtype1, module_defs[mtype0],
                     module_defs[mtype1])
            # append to "row" lists
            LUTs0.append(LUT)
            module_pair_efficiencies_vectors0.append(
                module_pair_efficiencies_vector)
        # append "rows" to double-level lists
        all_LUTs.append(LUTs0)
        all_module_pair_efficiencies_vectors.append(
            module_pair_efficiencies_vectors0)

    return petsird.DetectionEfficiencies(
        calibration_factor=calibration_factor,
        detection_bin_efficiencies=all_detection_bin_efficiencies,
        module_pair_sgidlut=all_LUTs,
        module_pair_efficiencies_vectors=all_module_pair_efficiencies_vectors,
    )


def get_scanner_info(
    module_defs: Iterable[CylindricalBlocksInfo]
) -> petsird.ScannerInformation:

    scanner_geometry = get_scanner_geometry(module_defs)

    # TODO scanner_info.bulk_materials

    # TOF info (in mm)
    allTofBinEdges = []
    allTofResolutions = []
    for mtype0, m_def0 in enumerate(module_defs):
        allTofBinEdgesThisModuleType = []
        allTofResolutionsThisModuleType = []
        for mtype1, m_def1 in enumerate(module_defs):
            # In this example, we determine the coincidence window from the radii.
            # Obviously, a real scanner likely will do something else.
            max_distance = m_def0.radius + m_def1.radius
            # somewhat arbitrary choice for TOF resolution between modules of different
            # types (but currently overwritten below anyway)
            tof_resolution_this_pair = (m_def0.tof_resolution +
                                        m_def1.tof_resolution) / 2
            if mtype0 == mtype1:
                num_tof_bins = m_def0.number_of_tof_bins
            else:
                # In this example, no TOF information between different module-types
                num_tof_bins = 1
                allTofResolutionsThisModuleType.append(1000)  # in mm
            tofBinEdgesThisPair = petsird.BinEdges(edges=numpy.linspace(
                -max_distance, max_distance, num_tof_bins +
                1, dtype="float32"))
            allTofBinEdgesThisModuleType.append(tofBinEdgesThisPair)
            allTofResolutionsThisModuleType.append(tof_resolution_this_pair)
        allTofBinEdges.append(allTofBinEdgesThisModuleType)
        allTofResolutions.append(allTofResolutionsThisModuleType)

    allEnergyBinEdges = []
    allEnergyResolutionsAt511 = []
    for m_def in module_defs:
        allEnergyBinEdges.append(
            petsird.BinEdges(
                edges=numpy.linspace(m_def.LLD,
                                     m_def.ULD,
                                     m_def.number_of_event_energy_bins + 1,
                                     dtype="float32")))
        allEnergyResolutionsAt511.append(
            m_def.energy_resolution)  # as fraction of 511

    # We need energy bin info before being able to construct the detection
    # efficiencies, so we first construct a scanner without the efficiencies
    scanner = petsird.ScannerInformation(
        model_name="PETSIRD_TEST",
        scanner_geometry=scanner_geometry,
        tof_bin_edges=allTofBinEdges,
        tof_resolution=allTofResolutions,
        event_energy_bin_edges=allEnergyBinEdges,
        energy_resolution_at_511=allEnergyResolutionsAt511)

    # Now add the efficiencies
    scanner.detection_efficiencies = get_detection_efficiencies(
        scanner, module_defs)

    scanner.coincidence_policy = petsird.CoincidencePolicy.REJECT_HIGHER_MULTIPLES
    scanner.single_events_are_stored = False
    scanner.prompt_events_are_stored = True
    scanner.delayed_events_are_stored = False
    scanner.triple_events_are_stored = False
    scanner.quadruple_events_are_stored = False

    return scanner


def get_header() -> petsird.Header:
    subject = petsird.Subject(id="123456")
    institution = petsird.Institution(
        name="Some institution",
        address="Didcot, Oxfordshire, OX11 0DE, UK",
    )
    return petsird.Header(
        exam=petsird.ExamInformation(subject=subject, institution=institution),
        scanner=get_scanner_info([mtype0_def, mtype1_def]),
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
            num_types_of_modules = header.scanner.scanner_geometry.number_of_module_types(
            )
            prompts_this_block = []
            for mtype0 in range(num_types_of_modules):
                prompts_mod0 = []
                for mtype1 in range(num_types_of_modules):
                    type_of_module_pair = petsird.TypeOfModulePair(
                        (mtype0, mtype1))
                    prompts_mod0.append(
                        list(
                            get_events(header, type_of_module_pair,
                                       num_prompts_this_block)))
                    prompts_this_block.append(prompts_mod0)
            # Normally we'd write multiple blocks, but here we have just one,
            # so let's write a tuple with just one element
            writer.write_time_blocks((petsird.TimeBlock.EventTimeBlock(
                petsird.EventTimeBlock(time_interval=time_interval,
                                       prompt_events=prompts_this_block)), ))
