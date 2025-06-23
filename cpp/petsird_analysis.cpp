/*
  Copyright (C) 2022-2023 Microsoft Corporation
  Copyright (C) 2023-2025 University College London

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
#include "petsird/helpers/geometry.h"
#include <xtensor/xview.hpp>
#include <xtensor/xio.hpp>
#include <iostream>
#include <variant>
#include <cstdlib>

void
print_usage_and_exit(char const* prog_name)
{
  std::cerr << "Usage:\n"
            << prog_name << " \\\n"
            << "    [--print_events] [--input petsird_filename]\n\n"
            << "Options:\n"
            << "    -e, --print-events: Print event info\n"
            << "    -i, --input       : Filename to read\n\n"
            << "Currently, the following (deprecated) usage is also allowed:\n"
            << prog_name << " [options] [---] petsird_filename\n"
            << "Use of '--' is then required if the filename starts with -\n";

  std::exit(EXIT_FAILURE);
}

int
main(int argc, char const* argv[])
{
  auto prog_name = argv[0];

  bool print_events = false;
  std::string filename;
  // option processing
  while (argc > 1 && (strncmp(argv[1], "-", 1) == 0))
    {
      if (strcmp(argv[1], "--print_events") == 0 || strcmp(argv[1], "-e") == 0)
        {
          print_events = true;
        }
      else if (strcmp(argv[1], "--input") == 0 || strcmp(argv[1], "-i") == 0)
        {
          filename = argv[2];
          ++argv;
          --argc;
        }
      else if (strcmp(argv[1], "--") == 0)
        {
          // next argument is filename
          ++argv;
          --argc;
          break;
        }
      else
        print_usage_and_exit(prog_name);
      ++argv;
      --argc;
    }
  // Check if the user has provided a file
  if (!filename.empty() && argc > 1)
    print_usage_and_exit(prog_name);
  else if (filename.empty())
    {
      if (argc < 2)
        print_usage_and_exit(prog_name);
      else
        filename = argv[1];
    }

  // Open the file
  PETSIRDReader reader(filename);
  petsird::Header header;
  reader.ReadHeader(header);

  std::cout << "Processing file: " << filename << std::endl;
  if (header.exam) // only do this if present
    std::cout << "Subject ID: " << header.exam->subject.id << std::endl;
  std::cout << "Types of modules: " << header.scanner.scanner_geometry.replicated_modules.size() << std::endl;
  // TODO more than 1 type of module
  const petsird::TypeOfModule type_of_module{ 0 };
  std::cout << "Number of modules of first type: "
            << header.scanner.scanner_geometry.replicated_modules[type_of_module].transforms.size() << std::endl;
  std::cout << "Number of elements in modules of first type: "
            << header.scanner.scanner_geometry.replicated_modules[type_of_module].object.detecting_elements.transforms.size()
            << std::endl;
  std::cout << "Total number of 'crystals' in modules of first type : "
            << petsird_helpers::get_num_det_els(header.scanner, type_of_module) << std::endl;
  // example printing of coordinates
  {
    const petsird::ExpandedDetectionBin expanded_detection_bin{ 0, 0, 0 };
    const auto box_shape = petsird_helpers::geometry::get_detecting_box(header.scanner, type_of_module, expanded_detection_bin);
    std::cout << "Coordinates of first detecting bin:\n";
    for (auto& corner : box_shape.corners)
      {
        std::cout << "  [" << corner.c[0] << ", " << corner.c[1] << ", " << corner.c[2] << "]\n";
      }
  }
  const auto& tof_bin_edges = header.scanner.tof_bin_edges[type_of_module][type_of_module];
  const auto num_tof_bins = tof_bin_edges.NumberOfBins();
  std::cout << "Number of TOF bins: " << num_tof_bins << std::endl;
  std::cout << "TOF bin edges: " << tof_bin_edges.edges << std::endl;
  const auto& event_energy_bin_edges = header.scanner.event_energy_bin_edges[type_of_module];
  const auto num_event_energy_bins = event_energy_bin_edges.NumberOfBins();
  std::cout << "Number of energy bins: " << num_event_energy_bins << std::endl;
  std::cout << "Event energy bin edges: " << event_energy_bin_edges.edges << std::endl;
  const auto energy_mid_points = (xt::view(event_energy_bin_edges.edges, xt::range(0, event_energy_bin_edges.edges.size() - 1))
                                  + xt::view(event_energy_bin_edges.edges, xt::range(1, event_energy_bin_edges.edges.size())))
                                 / 2;
  std::cout << "Event energy mid points: " << energy_mid_points << std::endl;

  std::cout << "Singles Histogram Level: ";
  switch (header.scanner.singles_histogram_level)
    {
    case petsird::SinglesHistogramLevelType::kNone:
      std::cout << "none\n";
      break;
    case petsird::SinglesHistogramLevelType::kModule:
      std::cout << "module\n";
      break;
    case petsird::SinglesHistogramLevelType::kAll:
      std::cout << "all\n";
      break;
    }
  if (header.scanner.singles_histogram_level != petsird::SinglesHistogramLevelType::kNone)
    {
      const auto& singles_histogram_energy_bin_edges = header.scanner.singles_histogram_energy_bin_edges[type_of_module];
      std::cout << "Singles Histogram Energy Bin Edges: " << singles_histogram_energy_bin_edges.edges << std::endl;
      std::cout << "Number of Singles Histogram Energy Windows: " << singles_histogram_energy_bin_edges.NumberOfBins()
                << std::endl;
    }

  petsird::TimeBlock time_block;

  // Process events in batches of up to 100
  float energy_1 = 0, energy_2 = 0;
  std::size_t num_prompts = 0;
  std::size_t num_delayeds = 0;
  float last_time = 0.F;
  while (reader.ReadTimeBlocks(time_block))
    {
      if (std::holds_alternative<petsird::EventTimeBlock>(time_block))
        {
          auto& event_time_block = std::get<petsird::EventTimeBlock>(time_block);
          last_time = event_time_block.time_interval.stop;
          const petsird::TypeOfModulePair type_of_module_pair{ 0, 0 };
          num_prompts += event_time_block.prompt_events[type_of_module_pair[0]][type_of_module_pair[1]].size();
          if (event_time_block.delayed_events)
            num_delayeds += (*event_time_block.delayed_events)[type_of_module_pair[0]][type_of_module_pair[1]].size();
          if (print_events)
            std::cout << "=====================  Prompt events in time block from " << last_time << " ==============\n";

          // Note: just doing one module-type ATM
          for (auto& event : event_time_block.prompt_events[type_of_module_pair[0]][type_of_module_pair[1]])
            {
              const auto expanded_detection_bin0
                  = petsird_helpers::expand_detection_bin(header.scanner, type_of_module_pair[0], event.detection_bins[0]);
              const auto expanded_detection_bin1
                  = petsird_helpers::expand_detection_bin(header.scanner, type_of_module_pair[1], event.detection_bins[1]);

              energy_1 += energy_mid_points[expanded_detection_bin0.energy_index];
              energy_2 += energy_mid_points[expanded_detection_bin1.energy_index];

              if (print_events)
                {
                  std::cout << "CoincidenceEvent(detectionBins=[" << event.detection_bins[0] << ", " << event.detection_bins[1]
                            << "], tofIdx=" << event.tof_idx << "])\n";
                  std::cout << "    "
                            << "[ExpandedDetectionBin(module=" << expanded_detection_bin0.module_index << ", "
                            << "el=" << expanded_detection_bin0.element_index << ", "
                            << "energy_index=" << expanded_detection_bin0.energy_index
                            << "), ExpandedDetectionBin(module=" << expanded_detection_bin1.module_index << ", "
                            << "el=" << expanded_detection_bin1.element_index << ", "
                            << "energy_index=" << expanded_detection_bin0.energy_index << ")]\n";
                  std::cout << "    efficiency:"
                            << petsird_helpers::get_detection_efficiency(header.scanner, type_of_module_pair, event) << "\n";
                }
            }
        }
    }

  std::cout << "Last time block at " << last_time << " ms\n";
  std::cout << "Number of prompt events: " << num_prompts << std::endl;
  std::cout << "Number of delayed events: " << num_delayeds << std::endl;
  if (num_prompts > 0)
    {
      std::cout << "Average energy_1: " << energy_1 / num_prompts << std::endl;
      std::cout << "Average energy_2: " << energy_2 / num_prompts << std::endl;
    }

  return 0;
}
