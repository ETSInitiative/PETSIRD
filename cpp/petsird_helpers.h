/*
  Copyright (C) 2024 University College London

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
      for (const auto& rep_volume : det_els)
        {
          num_det_els += rep_volume.transforms.size() * rep_module.transforms.size();
        }
    }
  return num_det_els;
}

struct ModuleAndElement
{
  int module;
  int el;
};

template <class T>
inline std::vector<ModuleAndElement>
get_module_and_element(const ScannerGeometry& scanner_geometry, const T& scanner_det_ids)
{
  assert(scanner_geometry.replicated_modules.size() == 1);
  const auto& rep_module = scanner_geometry.replicated_modules[0];
  assert(rep_module.object.detecting_elements.size() == 1);
  int num_modules = rep_module.transforms.size();

  std::vector<ModuleAndElement> result;
  for (int det : scanner_det_ids)
    {
      result.push_back({ det % num_modules, det / num_modules });
    }
  return result;
}

inline float
get_detection_efficiency(const ScannerInformation& scanner, const CoincidenceEvent& event)
{
  float eff = 1.0F;
  const auto& det_el_efficiencies = scanner.detection_efficiencies.det_el_efficiencies;
  if (!det_el_efficiencies)
    {
      eff *= ((*det_el_efficiencies)(event.detector_ids[0], event.energy_indices[0])
              * (*det_el_efficiencies)(event.detector_ids[1], event.energy_indices[1]));
    }
  const auto& module_pair_efficiencies_vector = scanner.detection_efficiencies.module_pair_efficiencies_vector;
  if (!module_pair_efficiencies_vector)
    {
      const auto& module_pair_SGID_LUT = scanner.detection_efficiencies.module_pair_sgidlut;
      assert(!module_pair_SGID_LUT);

      const auto mod_and_els = get_module_and_element(scanner.scanner_geometry, event.detector_ids);
      assert(scanner.scanner_geometry.replicated_modules.size() == 1);
      const int SGID = (*module_pair_SGID_LUT)(mod_and_els[0].module, mod_and_els[1].module);
      assert(SGID >= 0);

      const auto& module_pair_efficiencies = (*module_pair_efficiencies_vector)[SGID];
      assert(module_pair_efficiencies.sgid == static_cast<unsigned>(SGID));
      eff *= module_pair_efficiencies.values(mod_and_els[0].el, event.energy_indices[0], mod_and_els[1].el,
                                             event.energy_indices[1]);
    }
  return eff;
}

} // namespace petsird_helpers
