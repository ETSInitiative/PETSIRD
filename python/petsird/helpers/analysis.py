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
        num_module_types = len(scanner.scanner_geometry.replicated_modules)
        print("Types of modules: ", num_module_types)
        type_of_module = 0
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
                last_time = time_block.value.time_interval.stop
                for mtype0 in range(num_module_types):
                    for mtype1 in range(num_module_types):
                        mtype_pair = petsird.TypeOfModulePair((mtype0, mtype1))

                        # count events
                        prompt_events = time_block.value.prompt_events[mtype0][
                            mtype1]
                        num_prompts += len(prompt_events)
                        if time_block.value.delayed_events is not None:
                            num_delayeds += len(
                                time_block.value.delayed_events[mtype0]
                                [mtype1])

                        if print_events:
                            print(
                                "---------------------------- prompts for modules : ",
                                mtype_pair)

                        for event in prompt_events:
                            expanded_det_bin0 = expand_detection_bin(
                                scanner, mtype0, event.detection_bins[0])
                            expanded_det_bin1 = expand_detection_bin(
                                scanner, mtype1, event.detection_bins[1])

                            # accumulate energies to print average below
                            energy_1 += energy_mid_points[
                                expanded_det_bin0.energy_index]
                            energy_2 += energy_mid_points[
                                expanded_det_bin1.energy_index]
                            if print_events:
                                print(event)
                                print(
                                    "   ",
                                    expanded_det_bin0,
                                    ", ",
                                    expanded_det_bin1,
                                )
                                eff = get_detection_efficiency(
                                    scanner, mtype_pair, event)
                                print("    efficiency:", eff)

                                # get detector-element coordinates (relative to gantry)
                                box_shape0 = petsird.helpers.geometry.get_detecting_box(
                                    scanner, mtype0, expanded_det_bin0)
                                box_shape1 = petsird.helpers.geometry.get_detecting_box(
                                    scanner, mtype0, expanded_det_bin1)

                                # print some info
                                # (but not complete box, as it's a bit overwhelming)
                                print(
                                    "    mean of detection box 0:",
                                    sum([
                                        corner.c
                                        for corner in box_shape0.corners
                                    ]) / len(box_shape0.corners))
                                print(
                                    "    mean of detection box 1:",
                                    sum([
                                        corner.c
                                        for corner in box_shape1.corners
                                    ]) / len(box_shape1.corners))

        print(f"Last time block at {last_time} ms")
        print(f"Number of prompt events: {num_prompts}")
        print(f"Number of delayed events: {num_delayeds}")
        if num_prompts > 0:
            print(f"Average energy_1: {energy_1 / num_prompts}")
            print(f"Average energy_2: {energy_2 / num_prompts}")
