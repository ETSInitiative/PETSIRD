"""
Validation functions and helpers related to scanner geometry.
"""
from collections.abc import Sequence

import numpy as np
import numpy.typing as npt

import petsird


def is_clockwise_2D(points: Sequence[tuple[float, float]]) -> bool:
    """
    Determine if a polygon's vertices are ordered clockwise.

    Args:
        points: List of (x, y) tuples representing polygon vertices.
                Must contain at least 3 points.

    Returns:
        True if vertices are clockwise, False if counterclockwise.

    Raises:
        ValueError: If fewer than 3 points are provided.
    """
    if len(points) < 3:
        raise ValueError("A polygon must have at least 3 vertices.")

    # Calculate signed area using the shoelace formula
    # see e.g. https://mathworld.wolfram.com/PolygonArea.html
    signed_area = 0.0
    n = len(points)
    for i in range(n):
        x1, y1 = points[i]
        x2, y2 = points[(i + 1) % n]  # wrap around to first point
        signed_area += (x2 - x1) * (y2 + y1)

    # Clockwise if signed_area > 0 (based on this variant of shoelace method)
    return signed_area > 0


def is_clockwise_3D(points: Sequence[npt.NDArray], direction1: npt.NDArray,
                    direction2: npt.NDArray) -> bool:
    """check if projections of points on plane are in clockwise order

    The place is specified by 2 direction vectors.
    """
    return is_clockwise_2D([(np.dot(p, direction1), np.dot(p, direction2))
                            for p in points])


def is_lex_less(a: npt.NDArray, b: npt.NDArray) -> bool:
    """check if a <= b using lexicographic order """
    return np.all(a <= b)


def validate_detecting_element(box: petsird.BoxShape) -> None:
    """check constraints on a detecting element"""
    corners = [coordinate.c for coordinate in box.corners]

    if not len(corners) == 8:
        raise ValueError("detecting_element should have 8 vertices")

    d0 = np.array([1., 0., 0.])
    d1 = np.array([0., 1., 0.])
    d2 = np.array([0., 0., 1.])

    mean_first_face = sum(corners[0:4]) / 4
    mean_second_face = sum(corners[4:8]) / 4
    box_direction = mean_second_face - mean_first_face
    box_direction /= np.linalg.norm(box_direction)

    # check direction is roughly along first axis
    if np.dot(d0, box_direction) < np.cos(np.pi / 5):
        raise ValueError("detecting_element is not oriented correctly")

    # (incomplete) check if (0,0,0) (i.e. average detection point) is inside
    # TODO add more checks
    if np.dot(mean_first_face, box_direction) >= 0 or np.dot(
            mean_second_face, box_direction) <= 0:
        raise ValueError(
            "detecting_element doesn't seem to be contain the average detection point"
        )

    # check clockwise first face
    if not is_clockwise_3D(corners[0:3], d1, d2) or not is_clockwise_3D(
            corners[1:4], d1, d2):
        raise ValueError("detecting_element first face is not clock-wise")
    if not is_clockwise_3D(corners[4:7], d1, d2) or not is_clockwise_3D(
            corners[5:], d1, d2):
        raise ValueError("detecting_element first face is not clock-wise")

    # check order
    if not all(is_lex_less(corners[0], c) for c in corners[1:3]):
        raise ValueError(
            "detecting_element first face first corner is not 'smallest'")
    if not all(is_lex_less(corners[4], c) for c in corners[5:8]):
        raise ValueError(
            "detecting_element second face first corner is not 'smallest'")
