#  Copyright (C) 2022-2023 Microsoft Corporation
#  Copyright (C) 2023-2024 University College London
#
#  SPDX-License-Identifier: Apache-2.0

import sys
import numpy
import petsird
from petsird_helpers import (
    get_num_det_els,
    get_module_and_element,
    get_detection_efficiency,
)

if __name__ == "__main__":
    with petsird.BinaryPETSIRDReader(sys.stdin.buffer) as reader:
        header = reader.read_header()
        if header.exam is not None:
            print(f"Subject ID: {header.exam.subject.id}")
        print(f"Scanner name: { header.scanner.model_name}")
        print(
            f"Types of modules: {len(header.scanner.scanner_geometry.replicated_modules)}"
        )
        print(
            f"Number of modules of first type: {len(header.scanner.scanner_geometry.replicated_modules[0].transforms)}"
        )
        print(
            f"Number of types of detecting elements in modules of first type: {len(header.scanner.scanner_geometry.replicated_modules[0].object.detecting_elements)}"
        )
        print(
            f"Number of elements of first type in modules of first type: {len(header.scanner.scanner_geometry.replicated_modules[0].object.detecting_elements[0].transforms)}"
        )
        print(
            f"Total number of 'crystals' {get_num_det_els(header.scanner.scanner_geometry)}"
        )
        print(f"Number of TOF bins: {header.scanner.number_of_tof_bins()}")
        print(f"Number of energy bins: {header.scanner.number_of_energy_bins()}")

        energy_bin_edges = header.scanner.energy_bin_edges
        print(f"Energy bin edges: {energy_bin_edges}")
        energy_mid_points = (energy_bin_edges[:-1] + energy_bin_edges[1:]) / 2
        print(f"Energy mid points: {energy_mid_points}")
        print("SGID LUT:\n", header.scanner.detection_efficiencies.module_pair_sgidlut)
        energy_1, energy_2 = 0.0, 0.0
        num_prompts = 0
        last_time = 0
        for time_block in reader.read_time_blocks():
            if isinstance(time_block, petsird.TimeBlock.EventTimeBlock):
                last_time = time_block.value.start  # all TimeBlock types have this field
                num_prompts += len(time_block.value.prompt_events)
                print("=====================  Events in time block from ", last_time, " ==============")
                for event in time_block.value.prompt_events:
                    energy_1 += energy_mid_points[event.energy_indices[0]]
                    energy_2 += energy_mid_points[event.energy_indices[1]]

                    print(event)
                    print(
                        "   ",
                        get_module_and_element(
                            header.scanner.scanner_geometry, event.detector_ids
                        ),
                    )
                    print(
                        "    efficiency:", get_detection_efficiency(header.scanner, event)
                    )

        print(f"Last time block at {last_time} ms")
        print(f"Number of prompt events: {num_prompts}")
        print(f"Average energy_1: {energy_1 / num_prompts}")
        print(f"Average energy_2: {energy_2 / num_prompts}")
