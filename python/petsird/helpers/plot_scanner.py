#  Copyright (C) 2024 University College London
#
#  SPDX-License-Identifier: Apache-2.0

# basic plotting of the scanner geometry
# preliminary code!
import sys

import matplotlib.pyplot as plt
import numpy
import petsird
import petsird.helpers.geometry
from mpl_toolkits.mplot3d.art3d import Poly3DCollection


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
        det_els = rep_module.object.detecting_elements
        for mod_transform in rep_module.transforms:
            for transform in det_els.transforms:
                draw_BoxShape(
                    ax,
                    petsird.helpers.geometry.transform_BoxShape(
                        petsird.helpers.geometry.mult_transforms(
                            [mod_transform, transform]),
                        det_els.object.shape,
                    ),
                )
    plt.show()
