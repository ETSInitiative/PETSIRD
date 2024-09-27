/*
  Copyright (C) 2022-2023 Microsoft Corporation
  Copyright (C) 2023 University College London

  SPDX-License-Identifier: Apache-2.0
*/

#include "generated/hdf5/protocols.h"
#include <xtensor/xview.hpp>
#include <xtensor/xio.hpp>
#include <iostream>

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
  prd::hdf5::PrdExperimentReader reader(argv[1]);
  prd::Header header;
  reader.ReadHeader(header);

  std::cout << "Processing file: " << argv[1] << std::endl;
  if (header.exam) // only do this if present
    std::cout << "Subject ID: " << header.exam->subject.id << std::endl;
  // TODO std::cout << "Number of detectors: " << header.scanner.NumberOfDetectors() << std::endl;
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

  prd::TimeBlock time_block;

  // Process events in batches of up to 100
  float energy_1 = 0, energy_2 = 0;
  std::size_t num_prompts = 0;
  float last_time = 0.F;
  while (reader.ReadTimeBlocks(time_block))
    {
      last_time = time_block.id * header.scanner.listmode_time_block_duration;
      num_prompts += time_block.prompt_events.size();

      for (auto& event : time_block.prompt_events)
        {
          energy_1 += energy_mid_points[event.energy_indices[0]];
          energy_2 += energy_mid_points[event.energy_indices[1]];
        }
    }

  std::cout << "Last time block at " << last_time << " ms\n";
  std::cout << "Number of prompt events: " << num_prompts << std::endl;
  std::cout << "Average energy_1: " << energy_1 / num_prompts << std::endl;
  std::cout << "Average energy_2: " << energy_2 / num_prompts << std::endl;

  return 0;
}
