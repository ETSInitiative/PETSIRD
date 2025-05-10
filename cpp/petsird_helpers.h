/*
  Copyright (C) 2024, 2025 University College London

  SPDX-License-Identifier: Apache-2.0
*/

#include "generated/types.h"

namespace petsird_helpers
{

using namespace petsird;

inline std::size_t
get_num_det_els(const ScannerGeometry& scanner_geometry)
{
  std::size_t num_det_els = 0;
  for (const auto& rep_module : scanner_geometry.replicated_modules)
    {
      const auto& det_els = rep_module.object.detecting_elements;
      num_det_els += det_els.transforms.size() * rep_module.transforms.size();
    }
  return num_det_els;
}

struct ExpandedDetectionBin
{
  uint32_t module;
  uint32_t el;
  uint32_t energy_idx;
};

template <class T>
inline std::vector<ExpandedDetectionBin>
expand_detection_bins(const ScannerGeometry& scanner_geometry, const T& list_of_detection_bins)
{
  assert(scanner_geometry.replicated_modules.size() == 1);
  const auto& rep_module = scanner_geometry.replicated_modules[0];

  // TODO currently det_el_idx are uint32_t, so use this type to avoid compiler warnings
  const uint32_t num_el_per_module = rep_module.object.detecting_elements.transforms.size();

  std::vector<ExpandedDetectionBin> result;
  for (auto bin : list_of_detection_bins)
    {
      const auto& det = bin.det_el_idx;
      result.push_back({ det / num_el_per_module, det % num_el_per_module, bin.energy_idx });
    }
  return result;
}

inline float
get_detection_efficiency(const ScannerInformation& scanner, const CoincidenceEvent& event)
{
  float eff = 1.0F;
  const auto& detection_bin_efficiencies = scanner.detection_efficiencies.detection_bin_efficiencies;
  if (detection_bin_efficiencies)
    {
      eff *= ((*detection_bin_efficiencies)(event.detection_bins[0].det_el_idx, event.detection_bins[0].energy_idx)
              * (*detection_bin_efficiencies)(event.detection_bins[1].det_el_idx, event.detection_bins[1].energy_idx));
    }
  const auto& module_pair_efficiencies_vector = scanner.detection_efficiencies.module_pair_efficiencies_vector;
  if (module_pair_efficiencies_vector)
    {
      const auto& module_pair_SGID_LUT = scanner.detection_efficiencies.module_pair_sgidlut;
      assert(module_pair_SGID_LUT);

      const auto expanded_det_bins = expand_detection_bins(scanner.scanner_geometry, event.detection_bins);
      assert(scanner.scanner_geometry.replicated_modules.size() == 1);
      const int SGID = (*module_pair_SGID_LUT)(expanded_det_bins[0].module, expanded_det_bins[1].module);
      if (SGID < 0)
        {
          return 0.;
        }

      const auto& module_pair_efficiencies = (*module_pair_efficiencies_vector)[SGID];
      assert(module_pair_efficiencies.sgid == static_cast<unsigned>(SGID));
      eff *= module_pair_efficiencies.values(expanded_det_bins[0].el, event.detection_bins[0].energy_idx, expanded_det_bins[1].el,
                                             event.detection_bins[1].energy_idx);
    }
  return eff;
}

} // namespace petsird_helpers
