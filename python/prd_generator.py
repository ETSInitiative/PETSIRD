import sys
import math
import random

import prd


def get_header():
    offsets = [0.0, 1.0, 2.0]
    angles = [2 * math.pi * i / 10 for i in range(10)]
    detectors = [
        prd.Detector(offset=offset, angle=angle, id=i)
        for i, (offset, angle) in enumerate(
            (offset, angle) for offset in offsets for angle in angles
        )
    ]

    return prd.Header(
        subject=prd.Subject(id="123456"),
        institution=prd.Institution(
            name="Diamond Light Source",
            address="Harwell Science and Innovation Campus, Didcot, Oxfordshire, OX11 0DE, UK",
        ),
        detectors=detectors,
    )


def get_events(header: prd.Header):
    detector_count = header.number_of_detectors()
    for _ in range(1000):
        yield prd.CoincidenceEvent(
            detector_1_id=random.randint(0, detector_count),
            detector_2_id=random.randint(0, detector_count),
            energy_1=random.random(),
            energy_2=random.random(),
            delta_t=random.random(),
        )


if __name__ == "__main__":
    with prd.BinaryPrdExperimentWriter(sys.stdout.buffer) as writer:
        header = get_header()
        writer.write_header(header)
        writer.write_events(get_events(header))
