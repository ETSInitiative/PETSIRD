"""
py.test tests for some PETSIRD helpers.
"""
#  Copyright (C) 2026 University College London
#
#  SPDX-License-Identifier: Apache-2.0

import numpy as np


def test_is_lex_less():
    from petsird.validate.geometry import is_lex_less
    assert is_lex_less(np.array((1, 2, 3)), np.array((1, 3, 4)))
    assert is_lex_less(np.array((1, 2, 3)), np.array((4, 3, 3)))
    assert is_lex_less(np.array((1, 2, 3)), np.array((2, 2, 3)))


def test_is_clockwise_2D():
    from petsird.validate.geometry import is_clockwise_2D
    assert is_clockwise_2D([(0, 0), (0, 1), (1, 2), (.5, .1)])
    assert is_clockwise_2D(
        [np.array(p) for p in [(0, 0), (0, 1), (1, 2), (.5, .1)]])
    assert not is_clockwise_2D([(0, 0), (0, 1), (-1, 2), (-.5, .1)])


def test_is_clockwise_3D():
    from petsird.validate.geometry import is_clockwise_3D
    assert is_clockwise_3D(
        [np.array(p + (4, )) for p in [(0, 0), (0, 1), (1, 2), (.5, .1)]],
        np.array((1, 0, 0)), np.array((0, 1, 0)))
    assert is_clockwise_3D([
        np.array((p[0], 3, p[1])) for p in [(0, 0), (0, 1), (1, 2), (.5, .1)]
    ], np.array((1, 0, 0)), np.array((0, 0, 1)))
