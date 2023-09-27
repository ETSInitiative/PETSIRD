import sys
import numpy
import prd


if __name__ == "__main__":
    with prd.BinaryPrdExperimentReader(sys.stdin.buffer) as reader:
        header = reader.read_header()
        print(f"Subject ID: {header.subject.id}")
        print(f"Number of detectors: {header.scanner.number_of_detectors()}")

        energy_bin_edges = header.scanner.energy_bin_edges
        print(f"Energy bin edges: {energy_bin_edges}")
        energy_mid_points = (energy_bin_edges[:-1] + energy_bin_edges[1:])/2
        print(f"Energy mid points: {energy_mid_points}")

        coincidence_timing_lut = reader.read_coincidence_timing_lut()

        energy_1, energy_2 = 0.0, 0.0
        num_events = 0
        for event in reader.read_events():
            energy_1 += energy_mid_points[event.energy_1_idx]
            energy_2 += energy_mid_points[event.energy_2_idx]
            num_events += 1

        print(f"Number of events: {num_events}")
        print(f"Average energy_1: {energy_1 / num_events}")
        print(f"Average energy_2: {energy_2 / num_events}")
