/*
  Copyright (C) 2022-2023 Microsoft Corporation
  Copyright (C) 2023-2025 University College London

  SPDX-License-Identifier: Apache-2.0
*/

// (un)comment if you want HDF5 or binary output
#define USE_HDF5

#ifdef USE_HDF5
#  include "petsird/hdf5/protocols.h"
using petsird::hdf5::PETSIRDReader;
#else
#  include "petsird/binary/protocols.h"
using petsird::binary::PETSIRDReader;
#endif
#include "petsird_helpers.h"
#include "petsird_helpers/geometry.h"
#include <xtensor/xview.hpp>
#include <xtensor/xio.hpp>
#include <iostream>
#include <variant>
#include <cstdlib>
#include <vector>

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

// helper function to compute the mean of the corners in a BoxShape
petsird::Coordinate
mean_position(const petsird::BoxShape& box_shape)
{
  petsird::Coordinate mean;
  mean.c = { 0, 0, 0 };
  for (auto& corner : box_shape.corners)
    {
      mean.c += corner.c;
    }
  mean.c /= box_shape.corners.size();
  return mean;
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
  const auto& scanner = header.scanner;

  std::cout << "Processing file: " << filename << std::endl;
  if (header.exam) // only do this if present
    std::cout << "Subject ID: " << header.exam->subject.id << std::endl;
  const auto num_module_types = scanner.scanner_geometry.replicated_modules.size();
  std::cout << "Types of modules: " << num_module_types << std::endl;
  std::vector<yardl::NDArray<float, 1>> all_energy_mid_points;
  for (petsird::TypeOfModule type_of_module = 0; type_of_module < num_module_types; ++type_of_module)
    {
      std::cout << "------ Module type " << type_of_module << std::endl;
      std::cout << "Number of modules of this type: "
                << scanner.scanner_geometry.replicated_modules[type_of_module].transforms.size() << std::endl;
      std::cout << "Number of elements in modules of this type: "
                << scanner.scanner_geometry.replicated_modules[type_of_module].object.detecting_elements.transforms.size()
                << std::endl;
      std::cout << "Total number of 'crystals' in modules of this type : "
                << petsird_helpers::get_num_det_els(scanner, type_of_module) << std::endl;

      const auto& tof_bin_edges = scanner.tof_bin_edges[type_of_module][type_of_module];
      const auto num_tof_bins = tof_bin_edges.NumberOfBins();
      std::cout << "Number of TOF bins: " << num_tof_bins << std::endl;
      std::cout << "TOF bin edges: " << tof_bin_edges.edges << std::endl;
      const auto& event_energy_bin_edges = scanner.event_energy_bin_edges[type_of_module];
      const auto num_event_energy_bins = event_energy_bin_edges.NumberOfBins();
      std::cout << "Number of energy bins: " << num_event_energy_bins << std::endl;
      std::cout << "Event energy bin edges: " << event_energy_bin_edges.edges << std::endl;
      const auto energy_mid_points
          = (xt::view(event_energy_bin_edges.edges, xt::range(0, event_energy_bin_edges.edges.size() - 1))
             + xt::view(event_energy_bin_edges.edges, xt::range(1, event_energy_bin_edges.edges.size())))
            / 2;
      std::cout << "Event energy mid points: " << energy_mid_points << std::endl;
      all_energy_mid_points.push_back(energy_mid_points);

      std::cout << "Calibration factor: " << scanner.detection_efficiencies.calibration_factor << std::endl;
      std::cout << "Singles Histogram Level: ";
      switch (scanner.singles_histogram_level)
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
      if (scanner.singles_histogram_level != petsird::SinglesHistogramLevelType::kNone)
        {
          const auto& singles_histogram_energy_bin_edges = scanner.singles_histogram_energy_bin_edges[type_of_module];
          std::cout << "Singles Histogram Energy Bin Edges: " << singles_histogram_energy_bin_edges.edges << std::endl;
          std::cout << "Number of Singles Histogram Energy Windows: " << singles_histogram_energy_bin_edges.NumberOfBins()
                    << std::endl;
        }
    } // loop over type_of_module

  std::cout << "------------------------- " << std::endl;

  // Now read events and print some things
  petsird::TimeBlock time_block;
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

          if (print_events)
            std::cout << "=====================  Prompt events in time block from " << last_time << " ==============\n";

          for (unsigned mtype0 = 0; mtype0 < num_module_types; ++mtype0)
            {
              const auto& energy_mid_points0 = all_energy_mid_points[mtype0];
              for (unsigned mtype1 = 0; mtype1 < num_module_types; ++mtype1)
                {
                  const petsird::TypeOfModulePair mtype_pair{ mtype0, mtype1 };
                  const auto& energy_mid_points1 = all_energy_mid_points[mtype1];

                  // This code would need work to be able to handle a list-mode file without prompts
                  const auto& prompt_events = event_time_block.prompt_events[mtype0][mtype1];

                  // count events
                  num_prompts += prompt_events.size();
                  if (scanner.delayed_events_are_stored)
                    {
                      num_delayeds += event_time_block.delayed_events[mtype0][mtype1].size();
                    }

                  if (print_events)
                    {
                      std::cout << "---------------------------- prompts for modules : "
                                << "[" << mtype_pair[0] << ", " << mtype_pair[1] << "]\n";
                    }

                  for (const auto& event : prompt_events)
                    {
                      const auto expanded_det_bin0
                          = petsird_helpers::expand_detection_bin(scanner, mtype0, event.detection_bins[0]);
                      const auto expanded_det_bin1
                          = petsird_helpers::expand_detection_bin(scanner, mtype1, event.detection_bins[1]);

                      // accumulate energies to print average below
                      energy_1 += energy_mid_points0[expanded_det_bin0.energy_index];
                      energy_2 += energy_mid_points1[expanded_det_bin1.energy_index];

                      if (print_events)
                        {
                          std::cout << "CoincidenceEvent(detectionBins=[" << event.detection_bins[0] << ", "
                                    << event.detection_bins[1] << "], tofIdx=" << event.tof_idx << "])\n";
                          std::cout << "    "
                                    << "[ExpandedDetectionBin(module=" << expanded_det_bin0.module_index << ", "
                                    << "el=" << expanded_det_bin0.element_index << ", "
                                    << "energy_index=" << expanded_det_bin0.energy_index
                                    << "), ExpandedDetectionBin(module=" << expanded_det_bin1.module_index << ", "
                                    << "el=" << expanded_det_bin1.element_index << ", "
                                    << "energy_index=" << expanded_det_bin1.energy_index << ")]\n";

                          const auto eff = petsird_helpers::get_detection_efficiency(scanner, mtype_pair, event);
                          std::cout << "    efficiency: " << eff << std::endl;

                          const auto box_shape0
                              = petsird_helpers::geometry::get_detecting_box(scanner, mtype0, expanded_det_bin0);
                          const auto box_shape1
                              = petsird_helpers::geometry::get_detecting_box(scanner, mtype0, expanded_det_bin1);
                          const auto mean0 = mean_position(box_shape0);
                          const auto mean1 = mean_position(box_shape1);

                          std::cout << "    mean of detection box 0: [" << mean0.c[0] << ", " << mean0.c[1] << ", " << mean0.c[2]
                                    << "]\n";
                          std::cout << "    mean of detection box 1: [" << mean1.c[0] << ", " << mean1.c[1] << ", " << mean1.c[2]
                                    << "]\n";
                        }
                    }
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
