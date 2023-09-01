import sys
import prd


if __name__ == "__main__":
    with prd.BinaryPrdExperimentReader(sys.stdin.buffer) as reader:
        header = reader.read_header()
        print(f"Subject ID: {header.subject.id}")
        print(f"Number of detectors: {header.number_of_detectors()}")

        energy_1, energy_2 = 0.0, 0.0
        num_events = 0
        for event in reader.read_events():
            energy_1 += event.energy_1
            energy_2 += event.energy_2
            num_events += 1

        print(f"Number of events: {num_events}")
        print(f"Average energy_1: {energy_1 / num_events}")
        print(f"Average energy_2: {energy_2 / num_events}")
