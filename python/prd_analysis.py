import sys
import numpy
import prd


if __name__ == "__main__":
    with prd.BinaryPrdExperimentReader(sys.stdin.buffer) as reader:
        header = reader.read_header()
        print(f"Subject ID: {header.exam.subject.id}")
        print(f"Number of detectors: {header.scanner.number_of_detectors()}")
        print(f"Number of TOF bins: {header.scanner.number_of_tof_bins()}")
        print(f"Number of TOF bins: {header.scanner.number_of_energy_bins()}")

        energy_bin_edges = header.scanner.energy_bin_edges
        print(f"Energy bin edges: {energy_bin_edges}")
        energy_mid_points = (energy_bin_edges[:-1] + energy_bin_edges[1:])/2
        print(f"Energy mid points: {energy_mid_points}")

        energy_1, energy_2 = 0.0, 0.0
        num_events = 0
        last_time = 0
        for time_block in reader.read_time_blocks():
            last_time = time_block.id * header.scanner.listmode_time_block_duration
            for event in time_block.prompt_events:
                energy_1 += energy_mid_points[event.energy_1_idx]
                energy_2 += energy_mid_points[event.energy_2_idx]
                num_events += 1

        print(f"Last time block at {last_time} ms")
        print(f"Number of events: {num_events}")
        print(f"Average energy_1: {energy_1 / num_events}")
        print(f"Average energy_2: {energy_2 / num_events}")
