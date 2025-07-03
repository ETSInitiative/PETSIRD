/*
  Copyright (C) 2025 University College London

  SPDX-License-Identifier: Apache-2.0
*/

#ifndef __petsird_helpers_create_h__
#define __petsird_helpers_create_h__

#include "generated/types.h"
#include <array>

namespace petsird_helpers
{
namespace create
{

//! Helper function to create a std::vector<T>
/*! This function is added to have a 1D analogue construct_2D_nested_vector() */
template <typename T>
std::vector<T>
construct_vector(std::size_t size)
{
  std::vector<T> v(size);
  return v;
}

//! Helper function to create a std::vector<std::vector<T>> as a 2D array of size (size0, size1)
template <typename T>
std::vector<std::vector<T>>
construct_2D_nested_vector(std::size_t size0, std::size_t size1)
{
  std::vector<std::vector<T>> v(size0);
  for (auto& one_of_them : v)
    {
      one_of_them.resize(size1);
    }
  return v;
}

//! Set  various structures to have the correct size for the given num_types_of_modules
/*!
  This will set scanner.tof_bin_edges, scanner.tof_resolution,
  scanner.event_energy_bin_edges, scanner.energy_resolution_at_511, and
  and (optionally) scanner.detection_efficiencies.detection_bin_efficiencies,
  scanner.detection_efficiencies.module_pair_sgidlut and
  scanner.detection_efficiencies.module_pair_efficiencies_vectors
  to (nested) vectors of the appropriate type and size.

  Elements will be constructed via the default constructors, so you will still have
  to fill in the actual values.
*/
inline void
initialize_scanner_information_dimensions(petsird::ScannerInformation& scanner, const std::size_t num_module_types,
                                          bool allocate_detection_bin_efficiencies, bool allocate_module_pair_efficiencies)
{
  scanner.tof_bin_edges = construct_2D_nested_vector<petsird::BinEdges>(num_module_types, num_module_types);
  scanner.tof_resolution = construct_2D_nested_vector<float>(num_module_types, num_module_types);
  scanner.event_energy_bin_edges = construct_vector<petsird::BinEdges>(num_module_types);
  scanner.energy_resolution_at_511 = construct_vector<float>(num_module_types);

  scanner.detection_efficiencies = petsird::DetectionEfficiencies();

  if (allocate_detection_bin_efficiencies)
    {
      scanner.detection_efficiencies.detection_bin_efficiencies
          = construct_vector<petsird::DetectionBinEfficiencies>(num_module_types);
    }
  if (allocate_module_pair_efficiencies)
    {
      scanner.detection_efficiencies.module_pair_sgidlut
          = construct_2D_nested_vector<petsird::ModulePairSGIDLUT>(num_module_types, num_module_types);
      scanner.detection_efficiencies.module_pair_efficiencies_vectors
          = construct_2D_nested_vector<petsird::ModulePairEfficienciesVector>(num_module_types, num_module_types);
    }
}

} // namespace create
} // namespace petsird_helpers
#endif
