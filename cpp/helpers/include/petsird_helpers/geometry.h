/*
  Copyright (C) 2024, 2025 University College London

  SPDX-License-Identifier: Apache-2.0
*/
#ifndef __petsird_helpers_geometry_h__
#define __petsird_helpers_geometry_h__

#ifndef XTENSOR_VERSION_MAJOR
#  include <xtensor/xtensor.hpp>
#endif

#if XTENSOR_VERSION_MAJOR == 0 && XTENSOR_VERSION_MINOR < 26
#  include <xtensor/xarray.hpp>
#  include <xtensor/xview.hpp>
#  include <xtensor/xio.hpp>
#else
#  include <xtensor/containers/xarray.hpp>
#  include <xtensor/views/xview.hpp>
#  include <xtensor/io/xio.hpp>
#endif
#include <xtensor-blas/xlinalg.hpp>
#include <vector>
#include "petsird/types.h"
#include "petsird_helpers.h"
#include <algorithm>

namespace petsird_helpers
{
namespace geometry
{

// convenience functions for algebraic operations on Coordinate
inline petsird::Coordinate
operator+(const petsird::Coordinate& c1, const petsird::Coordinate& c2)
{
  petsird::Coordinate result;
  std::transform(c1.c.begin(), c1.c.end(), c2.c.begin(), result.c.begin(), std::plus<>{});
  return result;
}
inline petsird::Coordinate
operator-(const petsird::Coordinate& c1, const petsird::Coordinate& c2)
{
  petsird::Coordinate result;
  std::transform(c1.c.begin(), c1.c.end(), c2.c.begin(), result.c.begin(), std::minus<>{});
  return result;
}

inline petsird::Coordinate
operator+(const petsird::Coordinate& c, float arg)
{
  petsird::Coordinate result;
  std::transform(c.c.begin(), c.c.end(), result.c.begin(), [&](float a) { return a + arg; });
  return result;
}
inline petsird::Coordinate
operator-(const petsird::Coordinate& c, float arg)
{
  petsird::Coordinate result;
  std::transform(c.c.begin(), c.c.end(), result.c.begin(), [&](float a) { return a - arg; });
  return result;
}
inline petsird::Coordinate
operator*(const petsird::Coordinate& c, float arg)
{
  petsird::Coordinate result;
  std::transform(c.c.begin(), c.c.end(), result.c.begin(), [&](float a) { return a * arg; });
  return result;
}
inline petsird::Coordinate
operator/(const petsird::Coordinate& c, float arg)
{
  petsird::Coordinate result;
  std::transform(c.c.begin(), c.c.end(), result.c.begin(), [&](float a) { return a / arg; });
  return result;
}
//! Convert RigidTransformation to 4x4 matrix
inline xt::xarray<float>
transform_to_mat44(const petsird::RigidTransformation& transform)
{
  xt::xarray<float> bottom_row = { { 0.0f, 0.0f, 0.0f, 1.0f } };
  return xt::concatenate(xt::xtuple(transform.matrix, bottom_row), 0);
}

//! Convert 4x4 matrix to RigidTransformation
inline petsird::RigidTransformation
mat44_to_transform(const xt::xarray<float>& mat)
{
  petsird::RigidTransformation transform;
  transform.matrix = xt::view(mat, xt::range(0, 3), xt::all());
  return transform;
}

//! Convert Coordinate to homogeneous vector
inline xt::xarray<float>
coordinate_to_homogeneous(const petsird::Coordinate& coord)
{
  return xt::concatenate(xt::xtuple(coord.c, xt::xarray<float>{ 1.0f }), 0);
}

//! Convert homogeneous vector to Coordinate
inline petsird::Coordinate
homogeneous_to_coordinate(const xt::xarray<float>& hom_coord)
{
  petsird::Coordinate coord;
  coord.c = xt::view(hom_coord, xt::range(0, 3));
  return coord;
}

//! Multiply a list of transformations
inline petsird::RigidTransformation
mult_transforms(const std::vector<petsird::RigidTransformation>& transforms)
{
  xt::xarray<float> mat = xt::eye<float>(4);
  for (auto it = transforms.rbegin(); it != transforms.rend(); ++it)
    {
      mat = xt::linalg::dot(transform_to_mat44(*it), mat);
    }
  return mat44_to_transform(mat);
}

//! Apply transformations to a coordinate
inline petsird::Coordinate
mult_transforms_coord(const std::vector<petsird::RigidTransformation>& transforms, const petsird::Coordinate& coord)
{
  xt::xarray<float> hom = xt::linalg::dot(transform_to_mat44(mult_transforms(transforms)), coordinate_to_homogeneous(coord));
  return homogeneous_to_coordinate(hom);
}

//! Transform a BoxShape
inline petsird::BoxShape
transform_BoxShape(const petsird::RigidTransformation& transform, const petsird::BoxShape& box_shape)
{
  petsird::BoxShape new_box;
  for (size_t i = 0; i < box_shape.corners.size(); ++i)
    {
      new_box.corners[i] = mult_transforms_coord({ transform }, box_shape.corners[i]);
    }
  return new_box;
}

//! find the BoxShape corresponding to a ExpandedDetectionBin
inline petsird::BoxShape
get_detecting_box(const petsird::ScannerInformation& scanner, const petsird::TypeOfModule& type_of_module,
                  const petsird::ExpandedDetectionBin& expanded_detection_bin)
{
  const auto& rep_module = scanner.scanner_geometry.replicated_modules[type_of_module];
  const auto& det_els = rep_module.object.detecting_elements;
  const auto& mod_transform = rep_module.transforms[expanded_detection_bin.module_index];
  const auto& transform = det_els.transforms[expanded_detection_bin.element_index];
  return transform_BoxShape(mult_transforms({ mod_transform, transform }), det_els.object.shape);
}

//! return a cuboid volume
/*! The box will be centred around \a center_in_mm and align with the coordinate axes.
 */
inline petsird::BoxSolidVolume
construct_cuboid_box(std::array<float, 3> sizes_in_mm, std::array<float, 3> center_in_mm, unsigned material_id)
{
  using petsird::Coordinate;
  // shift
  const Coordinate shift{ { center_in_mm[0] - sizes_in_mm[0] / 2.F, center_in_mm[1] - sizes_in_mm[1] / 2.F,
                            center_in_mm[2] - sizes_in_mm[2] / 2.F } };
  petsird::BoxShape crystal_shape{ shift + Coordinate{ { 0, 0, 0 } },
                                   shift + Coordinate{ { 0, 0, sizes_in_mm[2] } },
                                   shift + Coordinate{ { 0, sizes_in_mm[1], sizes_in_mm[2] } },
                                   shift + Coordinate{ { 0, sizes_in_mm[1], 0 } },
                                   shift + Coordinate{ { sizes_in_mm[0], 0, 0 } },
                                   shift + Coordinate{ { sizes_in_mm[0], 0, sizes_in_mm[2] } },
                                   shift + Coordinate{ { sizes_in_mm[0], sizes_in_mm[1], sizes_in_mm[2] } },
                                   shift + Coordinate{ { sizes_in_mm[0], sizes_in_mm[1], 0 } } };

  petsird::BoxSolidVolume crystal{ crystal_shape, material_id };
  return crystal;
}

//! return a module of num_repetitions cuboids.
/*!
 Each detecting element will be a cuboid with the given \a sizes_in_mm.
 The cuboids will be spaced according to \a translations_in_mm.
 This allows having gaps between the cuboids.

 The average DOI for photons travelling along 1,0,0 of the first layer in the module
 will be at 0,0,0. This means that the module will be centered around y=0, z=0.

 Constraint: \a translations_in_mm >= sizes_in_mm.
*/
inline petsird::DetectorModule
construct_cuboid_module(std::array<float, 3> num_repetitions, std::array<float, 3> translations_in_mm,
                        std::array<float, 3> sizes_in_mm, float average_DOI, unsigned material_id = 1U)
{
  assert(translations_in_mm >= sizes_in_mm);

  petsird::ReplicatedBoxSolidVolume rep_volume;
  {
    // We need 0,0,0 of the box to correspond to the average_DOI position for
    // photons travelling along 1,0,0
    // This gives the location of the "incoming" plane:
    // -sizes_in_mm[0]/2 + center_in_mm[0] == -average_DOI
    rep_volume.object = construct_cuboid_box(sizes_in_mm, { sizes_in_mm[0] / 2.F - average_DOI, 0.F, 0.F }, material_id);
    const auto N0 = num_repetitions[0];
    const auto N1 = num_repetitions[1];
    const auto N2 = num_repetitions[2];
    for (int rep0 = 0; rep0 < N0; ++rep0)
      for (int rep1 = 0; rep1 < N1; ++rep1)
        for (int rep2 = 0; rep2 < N2; ++rep2)
          {
            petsird::RigidTransformation transform{ { { 1.F, 0.F, 0.F, rep0 * translations_in_mm[0] },
                                                      { 0.F, 1.F, 0.F, (rep1 - N1 / 2.F) * translations_in_mm[1] },
                                                      { 0.F, 0.F, 1.F, (rep2 - N2 / 2.F) * translations_in_mm[2] } } };
            rep_volume.transforms.push_back(transform);
          }
  }

  petsird::DetectorModule detector_module;
  detector_module.detecting_elements = std::move(rep_volume);

  return detector_module;
}

//! find the BoxShape corresponding to a DetectionBin
inline petsird::BoxShape
get_detecting_box(const petsird::ScannerInformation& scanner, const petsird::TypeOfModule& type_of_module,
                  const petsird::DetectionBin& detection_bin)
{
  return get_detecting_box(scanner, type_of_module, expand_detection_bin(scanner, type_of_module, detection_bin));
}

} // namespace geometry
} // namespace petsird_helpers
#endif
