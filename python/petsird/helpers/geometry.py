"""
Preliminary helpers for geometric calculations for PETSIRD data
"""

#  Copyright (C) 2024 - 2025 University College London
#
#  SPDX-License-Identifier: Apache-2.0

import numpy
import numpy.typing as npt
import petsird


def transform_to_mat44(
    transform: petsird.RigidTransformation, ) -> npt.NDArray[numpy.float32]:
    return numpy.vstack([transform.matrix, [0, 0, 0, 1]])


def mat44_to_transform(
        mat: npt.NDArray[numpy.float32]) -> petsird.RigidTransformation:
    return petsird.RigidTransformation(matrix=mat[0:3, :])


def coordinate_to_homogeneous(
        coord: petsird.Coordinate) -> npt.NDArray[numpy.float32]:
    return numpy.hstack([coord.c, 1])


def homogeneous_to_coordinate(
    hom_coord: npt.NDArray[numpy.float32], ) -> petsird.Coordinate:
    return petsird.Coordinate(c=hom_coord[0:3])


def mult_transforms(
    transforms: list[petsird.RigidTransformation],
) -> petsird.RigidTransformation:
    """multiply rigid transformations"""
    mat = numpy.array(
        ((1, 0, 0, 0), (0, 1, 0, 0), (0, 0, 1, 0), (0, 0, 0, 1)),
        dtype="float32",
    )

    for t in reversed(transforms):
        mat = numpy.matmul(transform_to_mat44(t), mat)
    return mat44_to_transform(mat)


def mult_transforms_coord(transforms: list[petsird.RigidTransformation],
                          coord: petsird.Coordinate) -> petsird.Coordinate:
    """apply list of transformations to coordinate"""
    # TODO better to multiply with coordinates in sequence, than first
    # multiplying the matrices
    hom = numpy.matmul(
        transform_to_mat44(mult_transforms(transforms)),
        coordinate_to_homogeneous(coord),
    )
    return homogeneous_to_coordinate(hom)


def transform_BoxShape(transform: petsird.RigidTransformation,
                       box_shape: petsird.BoxShape) -> petsird.BoxShape:
    return petsird.BoxShape(corners=[
        mult_transforms_coord([transform], c) for c in box_shape.corners
    ])
