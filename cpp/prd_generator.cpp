/*
  Copyright (C) 2022-2023 Microsoft Corporation
  Copyright (C) 2023-2024 University College London

  SPDX-License-Identifier: Apache-2.0
*/

#include <iostream>
#include <cmath>
#include <random>
#include "generated/hdf5/protocols.h"
#include <xtensor/xindex_view.hpp>

// these are constants for now
const uint32_t NUMBER_OF_ENERGY_BINS = 3;
const uint32_t NUMBER_OF_TOF_BINS = 300;
const float RADIUS = 400.F;
const std::array<float, 3> CRYSTAL_LENGTH{ 4.F, 4.F, 20.F };
const uint32_t NUMBER_OF_TIME_BLOCKS = 6;
const float COUNT_RATE = 500.F;

// single ring as example
prd::ScannerInformation
get_scanner_info()
{
  prd::DetectorModule detector_module;
  {
    // Define a cuboid
    using prd::Coordinate;
    prd::BoxShape crystal_shape{ Coordinate{ 0, 0, 0 },
                                 Coordinate{ 0, 0, CRYSTAL_LENGTH[2] },
                                 Coordinate{ 0, CRYSTAL_LENGTH[1], CRYSTAL_LENGTH[2] },
                                 Coordinate{ 0, CRYSTAL_LENGTH[1], 0 },
                                 Coordinate{ CRYSTAL_LENGTH[0], 0, 0 },
                                 Coordinate{ CRYSTAL_LENGTH[0], 0, CRYSTAL_LENGTH[2] },
                                 Coordinate{ CRYSTAL_LENGTH[0], CRYSTAL_LENGTH[1], CRYSTAL_LENGTH[2] },
                                 Coordinate{ CRYSTAL_LENGTH[0], CRYSTAL_LENGTH[1], 0 } };

    prd::SolidVolume crystal{ crystal_shape, /* material_id */ 1 };

    // Define a module of 1x2 cuboids
    prd::ReplicatedSolidVolume rep_volume;
    {
      rep_volume.solid_volume = crystal;
      // translate along first axis
      prd::RigidTransformation transform{ { { 1.F, 0.F, 0.F, RADIUS }, { 0.F, 1.F, 0.F, 0.F }, { 0.F, 0.F, 1.F, 0.F } } };
      rep_volume.transforms.push_back(transform);
      rep_volume.ids.push_back(0);
      // and along second axis
      transform.matrix(1,3) = CRYSTAL_LENGTH[1];
      rep_volume.transforms.push_back(transform);
      rep_volume.ids.push_back(1);
    }

    detector_module.detecting_elements.push_back(rep_volume);
    detector_module.detecting_element_ids.push_back(0);

  } // end detector_module

  prd::ScannerInformation scanner_info;
  {
    // fill in geometry
    {
      prd::ReplicatedDetectorModule rep_module;
      {
        // TODO why module_field?
        rep_module.module_field = detector_module;
        int module_id = 0;
        std::vector<float> angles;
        for (int i = 0; i < 10; ++i)
          {
            angles.push_back(static_cast<float>(2 * M_PI * i / 10));
          }
        for (auto angle : angles)
          {
            prd::RigidTransformation transform{ { { std::cos(angle), std::sin(angle), 0.F, 0.F },
                                                  { -std::sin(angle), std::cos(angle), 0.F, 0.F },
                                                  { 0.F, 0.F, 1.F, 0.F } } };
            rep_module.ids.push_back(module_id++);
            rep_module.transforms.push_back(transform);
          }
      }

      scanner_info.model_name = "PETSIRD_TEST";

      scanner_info.scanner_geometry.replicated_modules.push_back(rep_module);
      scanner_info.scanner_geometry.ids.push_back(0);
      // TODO scanner_info.bulk_materials

      // TOF and energy information
      {
        typedef yardl::NDArray<float, 1> FArray1D;
        // TOF info (in mm)
        FArray1D::shape_type tof_bin_edges_shape = { NUMBER_OF_TOF_BINS + 1 };
        FArray1D tof_bin_edges(tof_bin_edges_shape);
        for (std::size_t i = 0; i < tof_bin_edges.size(); ++i)
          tof_bin_edges[i] = (i - NUMBER_OF_TOF_BINS / 2.F) / NUMBER_OF_TOF_BINS * 2 * RADIUS;
        FArray1D::shape_type energy_bin_edges_shape = { NUMBER_OF_ENERGY_BINS + 1 };
        FArray1D energy_bin_edges(energy_bin_edges_shape);
        for (std::size_t i = 0; i < energy_bin_edges.size(); ++i)
          energy_bin_edges[i] = 430.F + i * (650.F - 430.F) / NUMBER_OF_ENERGY_BINS;
        scanner_info.tof_bin_edges = tof_bin_edges;
        scanner_info.tof_resolution = 9.4F; // in mm
        scanner_info.energy_bin_edges = energy_bin_edges;
        scanner_info.energy_resolution_at_511 = .11F;    // as fraction of 511
        scanner_info.listmode_time_block_duration = 1.F; // ms
      }
    }
  }
  return scanner_info;
}

prd::Header
get_header()
{
  prd::Subject subject;
  subject.id = "123456";
  prd::Institution institution;
  institution.name = "Diamond Light Source";
  institution.address = "Harwell Science and Innovation Campus, Didcot, Oxfordshire, OX11 0DE, UK";
  prd::ExamInformation exam_info;
  exam_info.subject = subject;
  exam_info.institution = institution;
  prd::Header header;
  header.exam = exam_info;
  header.scanner = get_scanner_info();
  return header;
}

// return pair of integers between 0 and max
std::pair<int, int>
get_random_pair(int max)
{
  int a = rand() % max;
  int b = rand() % max;
  return std::make_pair(a, b);
}

uint32_t
get_random_energy_value()
{
  return rand() % NUMBER_OF_ENERGY_BINS;
}

uint32_t
get_random_tof_value()
{
  return rand() % NUMBER_OF_TOF_BINS;
}

std::vector<prd::CoincidenceEvent>
get_events(const prd::Header&, std::size_t num_events)
{
  std::vector<prd::CoincidenceEvent> events;
  events.reserve(num_events);
  for (std::size_t i = 0; i < num_events; ++i)
    {
      const auto detectors = get_random_pair(1); // TODO header.scanner.NumberOfDetectors());
      prd::CoincidenceEvent e;
      e.detector_ids[0] = detectors.first;
      e.detector_ids[1] = detectors.second;
      e.energy_indices[0] = get_random_energy_value();
      e.energy_indices[1] = get_random_energy_value();
      e.tof_idx = get_random_tof_value();
      events.push_back(e);
    }
  return events;
}

int
main(int argc, char* argv[])
{
  // Check if the user has provided a file
  if (argc < 2)
    {
      std::cerr << "Please provide a filename to write to" << std::endl;
      return 1;
    }

  std::string outfile = argv[1];
  std::remove(outfile.c_str());
  prd::hdf5::PrdExperimentWriter writer(outfile);

  const auto header = get_header();
  writer.WriteHeader(header);

  std::random_device rd;
  std::mt19937 gen(rd());
  for (std::size_t t = 0; t < NUMBER_OF_TIME_BLOCKS; ++t)
    {
      std::poisson_distribution<> poisson(COUNT_RATE);
      const auto num_prompts_this_block = poisson(gen);
      const auto prompts_this_block = get_events(header, num_prompts_this_block);
      prd::TimeBlock time_block;
      time_block.id = t;
      time_block.prompt_events = prompts_this_block;
      writer.WriteTimeBlocks(time_block);
    }
  writer.EndTimeBlocks();

  // Check that we have completed protocol
  writer.Close();
  return 0;
}
