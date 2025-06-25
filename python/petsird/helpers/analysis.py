#  Copyright (C) 2022-2023 Microsoft Corporation
#  Copyright (C) 2023-2025 University College London
#
#  SPDX-License-Identifier: Apache-2.0

import argparse
import sys

import petsird
import petsird.helpers.geometry
from petsird.helpers import (expand_detection_bin, get_detection_efficiency,
                             get_num_det_els)


def parserCreator():
    parser = argparse.ArgumentParser(
        prog='petsird_analysis',
        description='Example program that lists basic content of a PETSIRD file'
    )
    parser.add_argument('-e', '--print_events', action='store_true')
    parser.add_argument(
        "-i",
        "--input",
        type=str,
        default=None,
        help="File to read from, or stdin if omitted",
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = parserCreator()
    file = None
    if args.input is None:
        file = sys.stdin.buffer
    else:
        file = open(args.input, "rb")
    print_events = args.print_events

    with petsird.BinaryPETSIRDReader(file) as reader:
        header = reader.read_header()
        scanner = header.scanner
        if header.exam is not None:
            print(f"Subject ID: {header.exam.subject.id}")
        print(f"Scanner name: {scanner.model_name}")
        type_of_module = 0
        print("Types of modules: ",
              len(scanner.scanner_geometry.replicated_modules))
        print(
            "Number of modules of first type: ",
            len(scanner.scanner_geometry.replicated_modules[type_of_module].
                transforms))
        print(
            "Number of elements in modules of first type: ",
            len(scanner.scanner_geometry.replicated_modules[type_of_module].
                object.detecting_elements.transforms))
        print("Total number of 'crystals': ",
              get_num_det_els(scanner, type_of_module))
        # example printing of coordinates
        expanded_detection_bin = petsird.ExpandedDetectionBin(module_index=0,
                                                              element_index=0,
                                                              energy_index=0)
        box_shape = petsird.helpers.geometry.get_detecting_box(
            header.scanner, type_of_module, expanded_detection_bin)
        print("Coordinates of first detecting bin:")
        for corner in box_shape.corners:
            print(corner.c)

        tof_bin_edges = scanner.tof_bin_edges[type_of_module][type_of_module]
        num_tof_bins = tof_bin_edges.number_of_bins()
        event_energy_bin_edges = scanner.event_energy_bin_edges[type_of_module]
        num_event_energy_bins = event_energy_bin_edges.number_of_bins()
        print("Number of TOF bins: ", num_tof_bins)
        print("Tof bin edges: ", tof_bin_edges)
        print("Number of energy bins: ", num_event_energy_bins)
        event_energy_bin_edges = event_energy_bin_edges.edges
        print("Event energy bin edges: ", event_energy_bin_edges)
        energy_mid_points = (event_energy_bin_edges[:-1] +
                             event_energy_bin_edges[1:]) / 2
        print("Event energy mid points: ", energy_mid_points)
        print("Singles histogram level: ", scanner.singles_histogram_level)
        if scanner.singles_histogram_level != petsird.SinglesHistogramLevelType.NONE:
            print(
                "Number of singles histograms energy windows for first module-type: ",
                scanner.singles_histogram_energy_bin_edges[type_of_module].
                number_of_bins())
            print("Singles histogram energy bin edges: ",
                  scanner.singles_histogram_energy_bin_edges)
        print("SGID LUT:\n",
              scanner.detection_efficiencies.module_pair_sgidlut)
        energy_1, energy_2 = 0.0, 0.0
        num_prompts = 0
        num_delayeds = 0
        last_time = 0
        for time_block in reader.read_time_blocks():
            if isinstance(time_block, petsird.TimeBlock.EventTimeBlock):
                # TODO just doing one module-type ATM
                type_of_module_pair = petsird.TypeOfModulePair((0, 0))
                last_time = time_block.value.time_interval.stop
                num_prompts += len(time_block.value.prompt_events[
                    type_of_module_pair[0]][type_of_module_pair[1]])
                if time_block.value.delayed_events is not None:
                    num_delayeds += len(time_block.value.delayed_events[
                        type_of_module_pair[0]][type_of_module_pair[1]])
                print("=====================  Events between module-types ",
                      type_of_module_pair, " in time block until ", last_time,
                      " ==============")
                for event in time_block.value.prompt_events[
                        type_of_module_pair[0]][type_of_module_pair[1]]:
                    expanded_detection_bin0 = expand_detection_bin(
                        scanner, type_of_module_pair[0],
                        event.detection_bins[0])
                    expanded_detection_bin1 = expand_detection_bin(
                        scanner, type_of_module_pair[1],
                        event.detection_bins[1])

                    energy_1 += energy_mid_points[
                        expanded_detection_bin0.energy_index]
                    energy_2 += energy_mid_points[
                        expanded_detection_bin1.energy_index]
                    if print_events:
                        print(event)
                        print(
                            "   ",
                            expanded_detection_bin0,
                            ", ",
                            expanded_detection_bin1,
                        )
                        print(
                            "    efficiency:",
                            get_detection_efficiency(scanner,
                                                     type_of_module_pair,
                                                     event))

        print(f"Last time block at {last_time} ms")
        print(f"Number of prompt events: {num_prompts}")
        print(f"Number of delayed events: {num_delayeds}")
        if num_prompts > 0:
            print(f"Average energy_1: {energy_1 / num_prompts}")
            print(f"Average energy_2: {energy_2 / num_prompts}")
