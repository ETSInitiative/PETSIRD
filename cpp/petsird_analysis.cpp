/*
  Copyright (C) 2022-2023 Microsoft Corporation
  Copyright (C) 2023-2024 University College London

  SPDX-License-Identifier: Apache-2.0
*/

// (un)comment if you want HDF5 or binary output
#define USE_HDF5

#ifdef USE_HDF5
#  include "generated/hdf5/protocols.h"
using petsird::hdf5::PETSIRDReader;
#else
#  include "generated/binary/protocols.h"
using petsird::binary::PETSIRDReader;
#endif
#include "petsird_helpers.h"
#include <xtensor/xview.hpp>
#include <xtensor/xio.hpp>
#include <iostream>
#include <variant>

int
main(int argc, char* argv[])
{
  // Check if the user has provided a file
  if (argc < 2)
    {
      std::cerr << "Please provide a file to read" << std::endl;
      return 1;
    }

  // Open the file
  PETSIRDReader reader(argv[1]);
  petsird::Header header;
  reader.ReadHeader(header);

  std::cout << "Processing file: " << argv[1] << std::endl;
  if (header.exam) // only do this if present
    std::cout << "Subject ID: " << header.exam->subject.id << std::endl;
  std::cout << "Types of modules: " << header.scanner.scanner_geometry.replicated_modules.size() << std::endl;
  std::cout << "Number of modules of first type: " << header.scanner.scanner_geometry.replicated_modules[0].transforms.size()
            << std::endl;
  std::cout << "Number of types of detecting elements in modules of first type: "
            << header.scanner.scanner_geometry.replicated_modules[0].object.detecting_elements.size() << std::endl;
  std::cout << "Number of elements of first type in modules of first type: "
            << header.scanner.scanner_geometry.replicated_modules[0].object.detecting_elements[0].transforms.size() << std::endl;
  std::cout << "Total number of 'crystals': " << petsird_helpers::get_num_det_els(header.scanner.scanner_geometry) << std::endl;

  std::cout << "Number of TOF bins: " << header.scanner.NumberOfTOFBins() << std::endl;
  std::cout << "Number of energy bins: " << header.scanner.NumberOfEnergyBins() << std::endl;

  const auto& tof_bin_edges = header.scanner.tof_bin_edges;
  std::cout << "TOF bin edges: " << tof_bin_edges << std::endl;
  const auto& energy_bin_edges = header.scanner.energy_bin_edges;
  std::cout << "Energy bin edges: " << energy_bin_edges << std::endl;
  const auto energy_mid_points = (xt::view(energy_bin_edges, xt::range(0, energy_bin_edges.size() - 1))
                                  + xt::view(energy_bin_edges, xt::range(1, energy_bin_edges.size())))
                                 / 2;
  std::cout << "Energy mid points: " << energy_mid_points << std::endl;

  petsird::TimeBlock time_block;

  // Process events in batches of up to 100
  float energy_1 = 0, energy_2 = 0;
  std::size_t num_prompts = 0;
  float last_time = 0.F;
  while (reader.ReadTimeBlocks(time_block))
    {
      if (std::holds_alternative<petsird::EventTimeBlock>(time_block))
        {
          auto& event_time_block = std::get<petsird::EventTimeBlock>(time_block);
          last_time = event_time_block.start;
          num_prompts += event_time_block.prompt_events.size();
          std::cout << "=====================  Events in time block from " << last_time << " ==============\n";

          for (auto& event : event_time_block.prompt_events)
            {
              energy_1 += energy_mid_points[event.energy_indices[0]];
              energy_2 += energy_mid_points[event.energy_indices[1]];

              std::cout << "CoincidenceEvent(detectorIds=[" << event.detector_ids[0] << ", " << event.detector_ids[1]
                        << "], tofIdx=" << event.tof_idx << ", energyIndices=[" << event.energy_indices[0] << ", "
                        << event.energy_indices[1] << "])\n";
              const auto module_and_elems
                  = petsird_helpers::get_module_and_element(header.scanner.scanner_geometry, event.detector_ids);
              std::cout << "    "
                        << "[ModuleAndElement(module=" << module_and_elems[0].module << ", "
                        << "el=" << module_and_elems[0].el << "), ModuleAndElement(module=" << module_and_elems[0].module << ", "
                        << "el=" << module_and_elems[0].el << ")]\n";
              std::cout << "    efficiency:" << petsird_helpers::get_detection_efficiency(header.scanner, event) << "\n";
            }
        }
    }

  std::cout << "Last time block at " << last_time << " ms\n";
  std::cout << "Number of prompt events: " << num_prompts << std::endl;
  std::cout << "Average energy_1: " << energy_1 / num_prompts << std::endl;
  std::cout << "Average energy_2: " << energy_2 / num_prompts << std::endl;

  return 0;
}
