/*
  Copyright (C) 2022-2023 Microsoft Corporation
  Copyright (C) 2023-2024 University College London

  SPDX-License-Identifier: Apache-2.0
*/

#include <iostream>
#include <cmath>
#include <random>
#include <algorithm>

// (un)comment if you want HDF5 or binary output
#define USE_HDF5

#ifdef USE_HDF5
#  include "generated/hdf5/protocols.h"
using petsird::hdf5::PETSIRDWriter;
#else
#  include "generated/binary/protocols.h"
using petsird::binary::PETSIRDWriter;
#endif

#include "petsird_helpers.h"

// these are constants for now
constexpr uint32_t NUMBER_OF_EVENT_ENERGY_BINS = 3;
constexpr uint32_t NUMBER_OF_TOF_BINS = 300;
constexpr float RADIUS = 400.F;
constexpr std::array<float, 3> CRYSTAL_LENGTH{ 20.F, 4.F, 4.F };
constexpr std::array<float, 3> NUM_CRYSTALS_PER_MODULE{ 2, 4, 7 };
constexpr uint32_t NUM_MODULES_ALONG_RING{ 20 };
constexpr uint32_t NUM_MODULES_ALONG_AXIS{ 2 };
constexpr float MODULE_AXIS_SPACING{ (NUM_CRYSTALS_PER_MODULE[2] + 4) * CRYSTAL_LENGTH[2] };

constexpr uint32_t NUMBER_OF_TIME_BLOCKS = 6;
constexpr float COUNT_RATE = 500.F;
constexpr float EVENT_TIME_BLOCK_DURATION = 1.F;

//! return a cuboid volume
petsird::BoxSolidVolume
get_crystal()
{
  using petsird::Coordinate;
  petsird::BoxShape crystal_shape{ Coordinate{ { 0, 0, 0 } },
                                   Coordinate{ { 0, 0, CRYSTAL_LENGTH[2] } },
                                   Coordinate{ { 0, CRYSTAL_LENGTH[1], CRYSTAL_LENGTH[2] } },
                                   Coordinate{ { 0, CRYSTAL_LENGTH[1], 0 } },
                                   Coordinate{ { CRYSTAL_LENGTH[0], 0, 0 } },
                                   Coordinate{ { CRYSTAL_LENGTH[0], 0, CRYSTAL_LENGTH[2] } },
                                   Coordinate{ { CRYSTAL_LENGTH[0], CRYSTAL_LENGTH[1], CRYSTAL_LENGTH[2] } },
                                   Coordinate{ { CRYSTAL_LENGTH[0], CRYSTAL_LENGTH[1], 0 } } };

  petsird::BoxSolidVolume crystal{ crystal_shape, /* material_id */ 1 };
  return crystal;
}

//! return a module of NUM_CRYSTALS_PER_MODULE cuboids
petsird::DetectorModule
get_detector_module()
{
  petsird::ReplicatedBoxSolidVolume rep_volume;
  {
    rep_volume.object = get_crystal();
    constexpr auto N0 = NUM_CRYSTALS_PER_MODULE[0];
    constexpr auto N1 = NUM_CRYSTALS_PER_MODULE[1];
    constexpr auto N2 = NUM_CRYSTALS_PER_MODULE[2];
    unsigned id = 0;
    for (int rep0 = 0; rep0 < N0; ++rep0)
      for (int rep1 = 0; rep1 < N1; ++rep1)
        for (int rep2 = 0; rep2 < N2; ++rep2)
          {
            petsird::RigidTransformation transform{ { { 1.0, 0.0, 0.0, RADIUS + rep0 * CRYSTAL_LENGTH[0] },
                                                      { 0.0, 1.0, 0.0, (rep1 - N1 / 2) * CRYSTAL_LENGTH[1] },
                                                      { 0.0, 0.0, 1.0, (rep2 - N2 / 2) * CRYSTAL_LENGTH[2] } } };
            rep_volume.transforms.push_back(transform);
            rep_volume.ids.push_back(id++);
          }
  }

  petsird::DetectorModule detector_module;
  detector_module.detecting_elements.push_back(rep_volume);
  detector_module.detecting_element_ids.push_back(0);

  return detector_module;
}

//! return scanner build by rotating a module around the (0,0,1) axis
petsird::ScannerGeometry
get_scanner_geometry()
{
  petsird::ReplicatedDetectorModule rep_module;
  {
    rep_module.object = get_detector_module();
    uint32_t module_id = 0;
    std::vector<float> angles;
    for (unsigned int i = 0; i < NUM_MODULES_ALONG_RING; ++i)
      {
        angles.push_back(static_cast<float>((2 * M_PI * i) / NUM_MODULES_ALONG_RING));
      }
    for (auto angle : angles)
      for (unsigned ax_mod = 0; ax_mod < NUM_MODULES_ALONG_AXIS; ++ax_mod)
        {
          petsird::RigidTransformation transform{ { { std::cos(angle), std::sin(angle), 0.F, 0.F },
                                                    { -std::sin(angle), std::cos(angle), 0.F, 0.F },
                                                    { 0.F, 0.F, 1.F, MODULE_AXIS_SPACING * ax_mod } } };
          rep_module.ids.push_back(module_id++);
          rep_module.transforms.push_back(transform);
        }
  }
  petsird::ScannerGeometry scanner_geometry;
  scanner_geometry.replicated_modules.push_back(rep_module);
  scanner_geometry.ids.push_back(0);
  return scanner_geometry;
}

petsird::DetectionEfficiencies
get_detection_efficiencies(const petsird::ScannerInformation& scanner)
{
  const auto num_det_els = petsird_helpers::get_num_det_els(scanner.scanner_geometry);
  petsird::DetectionEfficiencies detection_efficiencies;

  detection_efficiencies.detection_bin_efficiencies = xt::ones<float>({ num_det_els, scanner.NumberOfEventEnergyBins() });

  // only 1 type of module in the current scanner
  assert(scanner.scanner_geometry.replicated_modules.size() == 1);
  const auto& rep_module = scanner.scanner_geometry.replicated_modules[0];
  const auto num_modules = rep_module.transforms.size();

  // We will only use rotational symmetries (no translation along the axis yet)
  // We also assume all module-pairs are in coincidence, except those with the same angle.
  // Writing a module number as (z-position, angle):
  //   eff((z1,a1), (z2, a2)) == eff((z1,0), (z2, abs(a2-a1)))
  // or in linear indices
  //   eff(z1 + NZ * a1, z2 + NZ * a2) == eff(z1, z2 + NZ * abs(a2 - a1))
  // (coincident) SGIDs need to start from 0, so ignoring self-coincident angles
  constexpr auto num_SGIDs = NUM_MODULES_ALONG_AXIS * NUM_MODULES_ALONG_AXIS * (NUM_MODULES_ALONG_RING - 1);
  // SGID = z1 + NZ * (z2 + NZ * abs(a2 - a1) - 1)
  constexpr auto NZ = NUM_MODULES_ALONG_AXIS;
  detection_efficiencies.module_pair_sgidlut = yardl::NDArray<int, 2>({ num_modules, num_modules });
  auto& module_pair_SGID_LUT = *detection_efficiencies.module_pair_sgidlut;
  for (unsigned int mod1 = 0; mod1 < num_modules; ++mod1)
    {
      for (unsigned int mod2 = 0; mod2 < num_modules; ++mod2)
        {
          const auto z1 = mod1 % NZ;
          const auto a1 = mod1 / NZ;
          const auto z2 = mod2 % NZ;
          const auto a2 = mod2 / NZ;
          if (a1 == a2)
            {
              module_pair_SGID_LUT(mod1, mod2) = -1;
            }
          else
            {
              module_pair_SGID_LUT(mod1, mod2) = z1 + NZ * (z2 + NZ * (std::abs(int(a2) - int(a1)) - 1));
            }
        }
    }
  // assert(module_pair_SGID_LUT).max() == num_SGIDs - 1);

  // assign an empty vector first, and reserve correct size
  detection_efficiencies.module_pair_efficiencies_vector = petsird::ModulePairEfficienciesVector();
  detection_efficiencies.module_pair_efficiencies_vector->reserve(num_SGIDs);

  assert(rep_module.object.detecting_elements.size() == 1);
  const auto& detecting_elements = rep_module.object.detecting_elements[0];
  const auto num_det_els_in_module = detecting_elements.transforms.size();
  for (unsigned int SGID = 0; SGID < num_SGIDs; ++SGID)
    {
      // extract first module_pair for this SGID. However, as this currently unused, it is commented out
      // const auto& module_pair = *std::find(module_pair_SGID_LUT.begin(), module_pair_SGID_LUT.end(), SGID);
      petsird::ModulePairEfficiencies module_pair_efficiencies;
      module_pair_efficiencies.values = yardl::NDArray<float, 4>(
          { num_det_els_in_module, scanner.NumberOfEventEnergyBins(), num_det_els_in_module, scanner.NumberOfEventEnergyBins() });
      // give some (non-physical) value
      module_pair_efficiencies.values.fill(SGID);
      module_pair_efficiencies.sgid = SGID;
      detection_efficiencies.module_pair_efficiencies_vector->emplace_back(module_pair_efficiencies);
      assert(detection_efficiencies.module_pair_efficiencies_vector->size() == unsigned(SGID + 1));
    }

  return detection_efficiencies;
}

petsird::ScannerInformation
get_scanner_info()
{
  petsird::ScannerInformation scanner_info;
  scanner_info.model_name = "PETSIRD_TEST";

  scanner_info.scanner_geometry = get_scanner_geometry();

  // TODO scanner_info.bulk_materials

  // TOF and energy information
  {
    typedef yardl::NDArray<float, 1> FArray1D;
    // TOF info (in mm)
    FArray1D::shape_type tof_bin_edges_shape = { NUMBER_OF_TOF_BINS + 1 };
    FArray1D tof_bin_edges(tof_bin_edges_shape);
    for (std::size_t i = 0; i < tof_bin_edges.size(); ++i)
      tof_bin_edges[i] = (i - NUMBER_OF_TOF_BINS / 2.F) / NUMBER_OF_TOF_BINS * 2 * RADIUS;
    FArray1D::shape_type event_energy_bin_edges_shape = { NUMBER_OF_EVENT_ENERGY_BINS + 1 };
    FArray1D event_energy_bin_edges(event_energy_bin_edges_shape);
    for (std::size_t i = 0; i < event_energy_bin_edges.size(); ++i)
      event_energy_bin_edges[i] = 430.F + i * (650.F - 430.F) / NUMBER_OF_EVENT_ENERGY_BINS;
    scanner_info.tof_bin_edges = tof_bin_edges;
    scanner_info.tof_resolution = 9.4F; // in mm
    scanner_info.event_energy_bin_edges = event_energy_bin_edges;
    scanner_info.energy_resolution_at_511 = .11F; // as fraction of 511
  }

  scanner_info.detection_efficiencies = get_detection_efficiencies(scanner_info);

  scanner_info.coincidence_policy = petsird::CoincidencePolicy::kRejectMultiples;
  scanner_info.delayed_coincidences_are_stored = false;
  scanner_info.triple_events_are_stored = false;

  return scanner_info;
}

petsird::Header
get_header()
{
  petsird::Subject subject;
  subject.id = "123456";
  petsird::Institution institution;
  institution.name = "Diamond Light Source";
  institution.address = "Harwell Science and Innovation Campus, Didcot, Oxfordshire, OX11 0DE, UK";
  petsird::ExamInformation exam_info;
  exam_info.subject = subject;
  exam_info.institution = institution;
  petsird::Header header;
  header.exam = exam_info;
  header.scanner = get_scanner_info();
  return header;
}

// return uint between 0 and max
uint32_t
get_random_uint(int max)
{
  return static_cast<unsigned>(rand() % max);
}

uint32_t
get_random_energy_value()
{
  return rand() % NUMBER_OF_EVENT_ENERGY_BINS;
}

uint32_t
get_random_tof_value()
{
  return rand() % NUMBER_OF_TOF_BINS;
}

std::vector<petsird::CoincidenceEvent>
get_events(const petsird::Header& header, std::size_t num_events)
{
  std::vector<petsird::CoincidenceEvent> events;
  events.reserve(num_events);
  const auto num_det_els = petsird_helpers::get_num_det_els(header.scanner.scanner_geometry);
  for (std::size_t i = 0; i < num_events; ++i)
    {
      petsird::CoincidenceEvent e;
      e.detection_bins[0].energy_idx = get_random_energy_value();
      e.detection_bins[1].energy_idx = get_random_energy_value();
      // Generate random det_el_idxs until detection effficiency is not zero
      while (true)
        {
          e.detection_bins[0].det_el_idx = get_random_uint(num_det_els);
          e.detection_bins[1].det_el_idx = get_random_uint(num_det_els);
          if (petsird_helpers::get_detection_efficiency(header.scanner, e) > 0)
            {
              // in coincidence, we can get out of the loop
              break;
            }
        }
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
  PETSIRDWriter writer(outfile);

  const auto header = get_header();
  writer.WriteHeader(header);

  std::random_device rd;
  std::mt19937 gen(rd());
  for (std::size_t t = 0; t < NUMBER_OF_TIME_BLOCKS; ++t)
    {
      constexpr auto average_num = EVENT_TIME_BLOCK_DURATION * COUNT_RATE;
      std::poisson_distribution<> poisson(average_num);
      const auto num_prompts_this_block = poisson(gen);
      const auto prompts_this_block = get_events(header, num_prompts_this_block);
      petsird::EventTimeBlock time_block;
      time_block.time_interval.start = t * EVENT_TIME_BLOCK_DURATION;
      time_block.time_interval.stop = (t + 1) * EVENT_TIME_BLOCK_DURATION;
      time_block.prompt_events = prompts_this_block;
      writer.WriteTimeBlocks(time_block);
    }
  writer.EndTimeBlocks();

  // Check that we have completed protocol
  writer.Close();
  return 0;
}
