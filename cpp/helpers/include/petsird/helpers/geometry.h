#include <xtensor/xarray.hpp>
#include <xtensor/xview.hpp>
#include <xtensor/xio.hpp>
#include <xtensor-blas/xlinalg.hpp>
#include <vector>
#include "generated/types.h"

namespace petsird_helpers
{
namespace geometry
{

//! Convert RigidTransformation to 4x4 matrix
inline xt::xarray<float>
transform_to_mat44(const RigidTransformation& transform)
{
  xt::xarray<float> bottom_row = { { 0.0f, 0.0f, 0.0f, 1.0f } };
  return xt::concatenate(xt::xtuple(transform.matrix, bottom_row), 0);
}

//! Convert 4x4 matrix to RigidTransformation
inline RigidTransformation
mat44_to_transform(const xt::xarray<float>& mat)
{
  RigidTransformation transform;
  transform.matrix = xt::view(mat, xt::range(0, 3), xt::all());
  return transform;
}

//! Convert Coordinate to homogeneous vector
inline xt::xarray<float>
coordinate_to_homogeneous(const Coordinate& coord)
{
  return xt::concatenate(xt::xtuple(coord.c, xt::xarray<float>{ 1.0f }), 0);
}

//! Convert homogeneous vector to Coordinate
inline Coordinate
homogeneous_to_coordinate(const xt::xarray<float>& hom_coord)
{
  Coordinate coord;
  coord.c = xt::view(hom_coord, xt::range(0, 3));
  return coord;
}

//! Multiply a list of transformations
inline RigidTransformation
mult_transforms(const std::vector<RigidTransformation>& transforms)
{
  xt::xarray<float> mat = xt::eye<float>(4);
  for (auto it = transforms.rbegin(); it != transforms.rend(); ++it)
    {
      mat = xt::linalg::dot(transform_to_mat44(*it), mat);
    }
  return mat44_to_transform(mat);
}

//! Apply transformations to a coordinate
inline Coordinate
mult_transforms_coord(const std::vector<RigidTransformation>& transforms, const Coordinate& coord)
{
  xt::xarray<float> hom = xt::linalg::dot(transform_to_mat44(mult_transforms(transforms)), coordinate_to_homogeneous(coord));
  return homogeneous_to_coordinate(hom);
}

//! Transform a BoxShape
inline BoxShape
transform_BoxShape(const RigidTransformation& transform, const BoxShape& box_shape)
{
  BoxShape new_box;
  for (size_t i = 0; i < box_shape.corners.size(); ++i)
    {
      new_box.corners[i] = mult_transforms_coord({ transform }, box_shape.corners[i]);
    }
  return new_box;
}

//! find the BoxShape corresponding to a ExpandedDetectionBin
inline BoxShape
get_detecting_box(const ScannerInformation& scanner, const TypeOfModule& type_of_module,
                  const ExpandedDetectionBin& expanded_detection_bin)
{
  const auto& rep_module = scanner.scanner_geometry.replicated_modules[type_of_module];
  const auto& det_els = rep_module.object.detecting_elements;
  const auto& mod_transform = rep_module.transforms[expanded_detection_bin.module_index];
  const auto& transform = det_els.transforms[expanded_detection_bin.element_index];
  return transform_BoxShape(mult_transforms({ mod_transform, transform }), det_els.object.shape);
}

//! find the BoxShape corresponding to a DetectionBin
inline BoxShape
get_detecting_box(const ScannerInformation& scanner, const TypeOfModule& type_of_module, const DetectionBin& detection_bin)
{
  return get_detecting_box(scanner, type_of_module, expand_detection_bin(scanner, type_of_module, detection_bin));
}

} // namespace geometry
} // namespace petsird_helpers
