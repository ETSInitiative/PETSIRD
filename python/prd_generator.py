import sys
import math
import random
import numpy

import prd

# these are constants for now
NUMBER_OF_ENERGY_BINS = 3
NUMBER_OF_TOF_BINS = 300
RADIUS = 400
NUMBER_OF_TIME_FRAMES = 6
NUMBER_OF_EVENTS = 1000

# single ring as example
def get_scanner_info():
    radius = RADIUS
    angles = [2 * math.pi * i / 10 for i in range(10)]
    detectors = [
        prd.Detector(x=x, y=y, z=0, id=i)
        for i, (x,y) in enumerate(
            (radius*math.sin(angle), radius*math.cos(angle)) for angle in angles
        )
    ]
    # TOF info (in mm)
    tofBinEdges = numpy.linspace(-RADIUS, RADIUS, NUMBER_OF_TOF_BINS+1, dtype='float32')
    energyBinEdges = numpy.linspace(430, 650, NUMBER_OF_ENERGY_BINS+1, dtype='float32')
    return prd.ScannerInformation(
        detectors=detectors,
        tof_bin_edges=tofBinEdges,
        tof_resolution=9.4, # in mm
        energy_bin_edges = energyBinEdges,
        energy_resolution_at_511 = .11, # as fraction of 511
    )

def get_header():
    subject=prd.Subject(id="123456")
    institution=prd.Institution(
        name="Diamond Light Source",
        address="Harwell Science and Innovation Campus, Didcot, Oxfordshire, OX11 0DE, UK",
        )
    return prd.Header(
        exam = prd.ExamInformation(subject=subject, institution=institution),
        scanner=get_scanner_info(),
    )

def get_coincidence_timing_LUT():
    # 1 ms timeframes
    time_frames = [ prd.TimeInterval(start=t, stop=t+1) for t in range(NUMBER_OF_TIME_FRAMES)]
    lut = numpy.linspace(0, NUMBER_OF_EVENTS, NUMBER_OF_TIME_FRAMES, dtype='uint64')
    return prd.CoincidenceTimingLUT(
        time_frames = prd.TimeFrameInformation(time_frames=time_frames),
        event_numbers_at_time_frame_start = lut,
        )

def get_events(header: prd.Header):
    detector_count = header.scanner.number_of_detectors()
    for _ in range(NUMBER_OF_EVENTS):
        yield prd.CoincidenceEvent(
            detector_1_id=random.randrange(0, detector_count),
            detector_2_id=random.randrange(0, detector_count),
            energy_1_idx=random.randrange(0, NUMBER_OF_ENERGY_BINS),
            energy_2_idx=random.randrange(0, NUMBER_OF_ENERGY_BINS),
            tof_idx=random.randrange(0, NUMBER_OF_TOF_BINS),
        )


if __name__ == "__main__":
    with prd.BinaryPrdExperimentWriter(sys.stdout.buffer) as writer:
    #with prd.NDJsonPrdExperimentWriter(sys.stdout) as writer:
        header = get_header()
        writer.write_header(header)
        writer.write_coincidence_timing_lut(get_coincidence_timing_LUT())
        writer.write_events(get_events(header))
