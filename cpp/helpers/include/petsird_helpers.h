/*
  Copyright (C) 2024, 2025 University College London

  SPDX-License-Identifier: Apache-2.0
*/
#ifndef __petsird_helpers_h__
#define __petsird_helpers_h__

#include "generated/types.h"
#include <array>

namespace petsird_helpers
{

using namespace petsird;

// TODO remove?
inline std::size_t
get_num_det_els(const ScannerInformation& scanner, const TypeOfModule& type_of_module)
{
  const auto& rep_module = scanner.scanner_geometry.replicated_modules[type_of_module];
  const auto& det_els = rep_module.object.detecting_elements;
  return det_els.transforms.size() * rep_module.transforms.size();
}

//! Compute total number of detecting bins in all modules of the given type
inline std::size_t
get_num_detection_bins(const ScannerInformation& scanner, const TypeOfModule& type_of_module)
{
  const auto& rep_module = scanner.scanner_geometry.replicated_modules[type_of_module];
  const auto& det_els = rep_module.object.detecting_elements;
  const auto& energy_bin_edges = scanner.event_energy_bin_edges[type_of_module];
  return det_els.transforms.size() * rep_module.transforms.size() * energy_bin_edges.NumberOfBins();
}

//! Create a vector of ExpandedDetectionBins from a list/vector of DetectionBins
template <class T>
inline std::vector<ExpandedDetectionBin>
expand_detection_bins(const ScannerInformation& scanner, const TypeOfModule& type_of_module, const T& list_of_detection_bins)
{
  assert(type_of_module < scanner.scanner_geometry.replicated_modules.size());
  const auto& rep_module = scanner.scanner_geometry.replicated_modules[type_of_module];
  const auto& energy_bin_edges = scanner.event_energy_bin_edges[type_of_module];
  // TODO currently indices are uint32_t, so use this type to avoid compiler warnings
  const uint32_t num_en = energy_bin_edges.NumberOfBins();
  const uint32_t num_el_per_module = rep_module.object.detecting_elements.transforms.size();

  std::vector<ExpandedDetectionBin> result;
  for (auto bin : list_of_detection_bins)
    {
      const auto det = bin / num_en;
      result.push_back({ det / num_el_per_module, det % num_el_per_module, bin % num_en });
    }
  return result;
}

//! Expand a DetectionBin into a ExpandedDetectionBin
inline ExpandedDetectionBin
expand_detection_bin(const ScannerInformation& scanner, const TypeOfModule& type_of_module, const DetectionBin& detection_bin)
{
  // TODO very inefficient, but avoid re-implementation of code above.
  std::vector<DetectionBin> bins{ detection_bin };
  const auto expanded_bins = expand_detection_bins(scanner, type_of_module, bins);
  return expanded_bins[0];
}

//! Create a vector of DetectionBins from a list/vector of ExpandedDetectionBins
template <class T>
inline std::vector<DetectionBin>
make_detection_bins(const ScannerInformation& scanner, const TypeOfModule& type_of_module,
                    const T& list_of_expanded_detection_bins)
{
  assert(type_of_module < scanner.scanner_geometry.replicated_modules.size());
  const auto& rep_module = scanner.scanner_geometry.replicated_modules[type_of_module];
  const auto& energy_bin_edges = scanner.event_energy_bin_edges[type_of_module];
  const auto num_en = energy_bin_edges.NumberOfBins();
  const auto num_el_per_module = rep_module.object.detecting_elements.transforms.size();

  std::vector<DetectionBin> result;
  // TODO call reserve()
  for (auto bin : list_of_expanded_detection_bins)
    {
      // need to do static_cast to avoid compiler warning on narrowing from std::size_t to uint32_t
      const auto detection_bin
          = static_cast<DetectionBin>(bin.energy_index + (bin.element_index + bin.module_index * num_el_per_module) * num_en);
      result.push_back({ detection_bin });
    }
  return result;
}

//! Create a DetectionBin from a ExpandedDetectionBin
inline DetectionBin
make_detection_bin(const ScannerInformation& scanner, const TypeOfModule& type_of_module,
                   const ExpandedDetectionBin& expanded_detection_bin)
{
  // TODO very inefficient, but avoid re-implementation of code above.
  std::vector<ExpandedDetectionBin> expanded_detection_bins{ expanded_detection_bin };
  const auto detection_bins = make_detection_bins(scanner, type_of_module, expanded_detection_bins);
  return detection_bins[0];
}

inline float
get_detection_efficiency(const ScannerInformation& scanner, const TypeOfModulePair& type_of_module_pair,
                         const DetectionBin& detection_bin_1, const DetectionBin& detection_bin_2)
{
  float eff = 1.0F;
  const auto& detection_bin_efficiencies = scanner.detection_efficiencies.detection_bin_efficiencies;
  if (detection_bin_efficiencies)
    {
      eff *= ((*detection_bin_efficiencies)[type_of_module_pair[0]](detection_bin_1)
              * (*detection_bin_efficiencies)[type_of_module_pair[1]](detection_bin_2));
      if (eff == 0.F)
        return 0.F;
    }
  const auto& module_pair_efficiencies_vectors = scanner.detection_efficiencies.module_pair_efficiencies_vectors;
  if (module_pair_efficiencies_vectors)
    {
      assert(scanner.detection_efficiencies.module_pair_sgidlut);
      const auto& module_pair_SGID_LUT
          = (*scanner.detection_efficiencies.module_pair_sgidlut)[type_of_module_pair[0]][type_of_module_pair[1]];

      const auto expanded_det_bin0 = expand_detection_bin(scanner, type_of_module_pair[0], detection_bin_1);
      const auto expanded_det_bin1 = expand_detection_bin(scanner, type_of_module_pair[1], detection_bin_2);
      const int SGID = module_pair_SGID_LUT(expanded_det_bin0.module_index, expanded_det_bin1.module_index);
      if (SGID < 0)
        {
          return 0.;
        }

      const auto& module_pair_efficiencies
          = (*module_pair_efficiencies_vectors)[type_of_module_pair[0]][type_of_module_pair[1]][SGID];
      assert(module_pair_efficiencies.sgid == static_cast<unsigned>(SGID));
      const auto& num_en0 = scanner.event_energy_bin_edges[type_of_module_pair[0]].NumberOfBins();
      const auto& num_en1 = scanner.event_energy_bin_edges[type_of_module_pair[1]].NumberOfBins();
      // TODO create helper for next calculation
      eff *= module_pair_efficiencies.values(expanded_det_bin0.element_index * num_en0 + expanded_det_bin0.energy_index,
                                             expanded_det_bin1.element_index * num_en1 + expanded_det_bin1.energy_index);
    }
  return eff;
}

inline float
get_detection_efficiency(const ScannerInformation& scanner, const TypeOfModulePair& type_of_module_pair,
                         const CoincidenceEvent& event)
{
  return get_detection_efficiency(scanner, type_of_module_pair, event.detection_bins[0], event.detection_bins[1]);
}

} // namespace petsird_helpers

#endif
