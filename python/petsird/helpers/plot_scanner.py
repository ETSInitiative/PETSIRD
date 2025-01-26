#  Copyright (C) 2024 University College London
#
#  SPDX-License-Identifier: Apache-2.0

# basic plotting of the scanner geometry
# preliminary code!
import sys

import matplotlib.pyplot as plt
import numpy
import numpy.typing as npt
import petsird
from mpl_toolkits.mplot3d.art3d import Poly3DCollection


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


def draw_BoxShape(ax, box: petsird.BoxShape) -> None:
    vertices = numpy.array([c.c for c in box.corners])
    edges = [
        [vertices[j] for j in [0, 1, 2, 3]],
        [vertices[j] for j in [4, 5, 6, 7]],
        [vertices[j] for j in [0, 1, 5, 4]],
        [vertices[j] for j in [2, 3, 7, 6]],
        [vertices[j] for j in [1, 2, 6, 5]],
        [vertices[j] for j in [4, 7, 3, 0]],
    ]
    box = Poly3DCollection(edges, alpha=0.25, linewidths=1, edgecolors="r")
    ax.add_collection3d(box)


if __name__ == "__main__":
    reader = petsird.BinaryPETSIRDReader(sys.stdin.buffer)
    header = reader.read_header()

    # Create a new figure
    fig = plt.figure()

    # Add a 3D subplot
    ax = fig.add_subplot(111, projection="3d")

    # draw all crystals
    for rep_module in header.scanner.scanner_geometry.replicated_modules:
        det_el = rep_module.object.detecting_elements
        for mod_transform in rep_module.transforms:
            for rep_volume in det_el:
                for transform in rep_volume.transforms:
                    draw_BoxShape(
                        ax,
                        transform_BoxShape(
                            mult_transforms([mod_transform, transform]),
                            rep_volume.object.shape,
                        ),
                    )
    plt.show()
