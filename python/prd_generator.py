import sys
import math
import random
import numpy
from collections.abc import Iterator
import prd

# these are constants for now
NUMBER_OF_ENERGY_BINS = 3
NUMBER_OF_TOF_BINS = 300
RADIUS = 400
CRYSTAL_LENGTH = (4, 4, 20)
NUMBER_OF_TIME_BLOCKS = 6
NUMBER_OF_EVENTS = 1000
COUNT_RATE = 500


# single ring as example
def get_scanner_info() -> prd.ScannerInformation:
    crystal_shape = prd.BoxShape(
        corners=[
            prd.Coordinate(c=(0, 0, 0)),
            prd.Coordinate(c=(0, 0, CRYSTAL_LENGTH[2])),
            prd.Coordinate(c=(0, CRYSTAL_LENGTH[1], CRYSTAL_LENGTH[2])),
            prd.Coordinate(c=(0, CRYSTAL_LENGTH[1], 0)),
            prd.Coordinate(c=(CRYSTAL_LENGTH[0], 0, 0)),
            prd.Coordinate(c=(CRYSTAL_LENGTH[0], 0, CRYSTAL_LENGTH[2])),
            prd.Coordinate(c=(CRYSTAL_LENGTH[0], CRYSTAL_LENGTH[1], CRYSTAL_LENGTH[2])),
            prd.Coordinate(c=(CRYSTAL_LENGTH[0], CRYSTAL_LENGTH[1], 0)),
        ]
    )
    # handle unions in Python
    crystal_shape = prd.GeometricShape(
        shape=prd.BoxShapeOrAnnulusShape.BoxShape(crystal_shape)
    )
    crystal = prd.SolidVolume(shape=crystal_shape, material_id=1)

    # Define a module of 1x2 cuboids
    rep_volume = prd.ReplicatedSolidVolume(solid_volume=crystal)
    # translate along first axis
    transform = prd.RigidTransformation(
        matrix=numpy.array(
            ((1.0, 0.0, 0.0, RADIUS), (0.0, 1.0, 0.0, 0.0), (0.0, 0.0, 1.0, 0.0)),
            dtype="float32",
        )
    )
    rep_volume.transforms.append(transform)
    rep_volume.ids.append(0)
    # and along second axis
    # TODO does not work: TypeError: 'RigidTransformation' object does not support item assignment
    # transform[1, 3] = CRYSTAL_LENGTH[1]
    rep_volume.transforms.append(transform)
    rep_volume.ids.append(1)

    detector_module = prd.DetectorModule(
        detecting_elements=[rep_volume], detecting_element_ids=[0]
    )

    radius = RADIUS
    angles = [2 * math.pi * i / 10 for i in range(10)]

    rep_module = prd.ReplicatedDetectorModule(module=detector_module)
    module_id = 0
    for angle in angles:
        transform = prd.RigidTransformation(
            matrix=numpy.array(
                (
                    (math.cos(angle), math.sin(angle), 0.0, 0.0),
                    (-math.sin(angle), math.cos(angle), 0.0, 0.0),
                    (0.0, 0.0, 1.0, 0.0),
                ),
                dtype="float32",
            )
        )
        rep_module.ids.append(module_id)
        module_id += 1
        rep_module.transforms.append(transform)

    scanner_geometry = prd.ScannerGeometry(replicated_modules=[rep_module], ids=[0])
    # TODO scanner_info.bulk_materials

    # TOF info (in mm)
    tofBinEdges = numpy.linspace(
        -RADIUS, RADIUS, NUMBER_OF_TOF_BINS + 1, dtype="float32"
    )
    energyBinEdges = numpy.linspace(
        430, 650, NUMBER_OF_ENERGY_BINS + 1, dtype="float32"
    )
    return prd.ScannerInformation(
        model_name="PETSIRD_TEST",
        scanner_geometry=scanner_geometry,
        tof_bin_edges=tofBinEdges,
        tof_resolution=9.4,  # in mm
        energy_bin_edges=energyBinEdges,
        energy_resolution_at_511=0.11,  # as fraction of 511
        listmode_time_block_duration=1,  # ms
    )


def get_header() -> prd.Header:
    subject = prd.Subject(id="123456")
    institution = prd.Institution(
        name="Diamond Light Source",
        address="Harwell Science and Innovation Campus, Didcot, Oxfordshire, OX11 0DE, UK",
    )
    return prd.Header(
        exam=prd.ExamInformation(subject=subject, institution=institution),
        scanner=get_scanner_info(),
    )


def get_events(header: prd.Header, num_events: int) -> Iterator[prd.CoincidenceEvent]:
    detector_count = 5  # TODO header.scanner.number_of_detectors()
    for _ in range(num_events):
        yield prd.CoincidenceEvent(
            detector_ids=[
                random.randrange(0, detector_count),
                random.randrange(0, detector_count),
            ],
            energy_indices=[
                random.randrange(0, NUMBER_OF_ENERGY_BINS),
                random.randrange(0, NUMBER_OF_ENERGY_BINS),
            ],
            tof_idx=random.randrange(0, NUMBER_OF_TOF_BINS),
        )


if __name__ == "__main__":
    # numpy random number generator
    rng = numpy.random.default_rng()

    with prd.BinaryPrdExperimentWriter(sys.stdout.buffer) as writer:
        # with prd.NDJsonPrdExperimentWriter(sys.stdout) as writer:
        header = get_header()
        writer.write_header(header)
        for t in range(NUMBER_OF_TIME_BLOCKS):
            num_prompts_this_block = rng.poisson(COUNT_RATE)
            prompts_this_block = list(get_events(header, num_prompts_this_block))
            # Normally we'd write multiple blocks, but here we have just one, so let's write a tuple with just one element
            writer.write_time_blocks(
                (prd.TimeBlock(id=t, prompt_events=prompts_this_block),)
            )
