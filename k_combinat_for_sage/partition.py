# -*- coding: utf-8 -*-
r"""
Sage has a builtin `Partition <https://doc.sagemath.org/html/en/reference/combinat/sage/combinat/partition.html>`_ object.  *This* module adds extra useful functions for partitions:

AUTHORS:

- Matthew Lancellotti (2018): Initial version
"""

#*****************************************************************************
#  Copyright (C) 2018 Matthew Lancellotti <mvlancellotti@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.all import *
# ^*^ sphinx insert ^*^

# HELPERS


def is_weakly_decreasing(li):
    r""" Return whether every term in the iterable ``li`` is greater than or equal to the following term.

    Used internally to check when a :class:`Partition`, :class:`Composition`, or list is weakly decreasing.

    EXAMPLES::

        sage: is_weakly_decreasing([3, 2, 1])
        True
        sage: is_weakly_decreasing([3, 2, 2])
        True
        sage: is_weakly_decreasing([3, 2, 3])
        False
    """
    return all(li[i] >= li[i+1] for i in range(len(li)-1))


def is_strictly_decreasing(li):
    r""" Return whether every term in the iterable ``li`` is greater than the following term.

    Used internally to check when a :class:`Partition`, :class:`Composition`, or list is strictly decreasing.

    EXAMPLES::

        sage: is_strictly_decreasing([3, 2, 1])
        True
        sage: is_strictly_decreasing([3, 2, 2])
        False
        sage: is_strictly_decreasing([3, 2, 3])
        False
    """
    return all(li[i] > li[i+1] for i in range(len(li)-1))


def _is_sequence(obj):
    r""" Helper function for internal use.

    Return whether ``obj`` is one of our allowed 'compositions'.

    EXAMPLES::

        sage: _is_sequence([3, 2, 2])
        True
        sage: _is_sequence(Composition([3, 2, 2]))
        True
        sage: _is_sequence(Partition([3, 2, 2]))
        True
        sage: _is_sequence(vector([3, 2, 2]))
        False
    """
    return isinstance(obj, (list, Composition, Partition))


def k_rectangle_dimension_list(k):
    r""" Return the list of dimension pairs `(h, w)` such that `h + w = k + 1`.

    This exists mainly as a helper function for :meth:`partition.has_rectangle` and :meth:`k_shape.is_reducible`.

    EXAMPLES::

        sage: k_rectangle_dimension_list(3)
        [(3, 1), (2, 2), (1, 3)]
    """
    return [(k-i+1, i) for i in range(1, k+1)]


# Partition stuff
def k_size(ptn, k):
    r""" Given a partition ``ptn`` and a ``k``, return the size of the `k`-boundary.

    This is the same as the length method :meth:`sage.combinat.core.Core.length` of the :class:`sage.combinat.core.Core` object, with the exception that here we don't require ``ptn`` to be a `k+1`-core.

    EXAMPLES::

        sage: Partition([2, 1, 1]).k_size(1)
        2
        sage: Partition([2, 1, 1]).k_size(2)
        3
        sage: Partition([2, 1, 1]).k_size(3)
        3
        sage: Partition([2, 1, 1]).k_size(4)
        4

    ..  SEEALSO::

        :meth:`k_boundary`, :meth:`SkewPartition.size`
    """
    ptn = Partition(ptn)
    return ptn.k_boundary(k).size()


def boundary(ptn):
    r""" Return the integer coordinates of points on the boundary of ``ptn``.

    The boundary of a partition is the set `\{ \text{NE}(d) \mid \forall d\:\text{diagonal} \}`.  That is, for every diagonal line `y = x + b` where `b \in \mathbb{Z}`, we find the northeasternmost (NE) point on that diagonal which is also in the Ferrer's diagram (here, the Ferrer's diagram is interpreted as 1 x 1 cells in the Euclidean plane).

    The boundary will go from bottom-right to top-left in the French convention.

    EXAMPLES:

    The partition (1)

    .. image:: _static/boundary-1.JPG
        :height: 140px
        :align: center
        :alt: The shape of partition 1 on a cartesian plane with the points (1, 0), (1, 1), and (0, 1) labelled.

    has boundary [(1, 0), (1, 1), (0, 1)]::

        sage: boundary(Partition([1]))
        [(1,0), (1,1), (0,1)]

    The partition (3, 1)

    .. image:: _static/boundary-2.JPG
        :height: 170px
        :align: center
        :alt: The shape of partition (3, 1) on a cartesian plane with the points (3, 0), (3, 1), (2, 1), (1, 1), (1, 2), and (0, 2) labelled.

    has boundary [(3, 0), (3, 1), (2, 1), (1, 1), (1, 2), (0, 2)]::

        sage: boundary(Partition([3, 1]))
        [(3,0), (3,1), (2,1), (1,1), (1,2), (0,2)]

    ..  SEEALSO::

        :meth:`k_rim`.  You might have been looking for :meth:`k_boundary` instead.
    """
    def horizontal_piece(xy, bdy):
        (start_x, start_y) = xy
        if not bdy:
            h_piece = [(start_x, start_y)]
        else:
            stop_x = bdy[-1][0]
            y = start_y  # y never changes
            h_piece = [(x, y) for x in range(start_x, stop_x)]
        h_piece = list(reversed(h_piece))
        return h_piece
    bdy = []
    for i, part in enumerate(ptn):
        (cell_x, cell_y) = (part - 1, i)
        (x, y) = (cell_x + 1, cell_y + 1)
        bdy += horizontal_piece((x, y - 1), bdy)
        bdy.append((x, y))
    # add final "top-left" horizontal piece
    (top_left_x, top_left_y) = (0, len(ptn))
    bdy += horizontal_piece((top_left_x, top_left_y), bdy)
    return bdy


def k_rim(ptn, k):
    r""" Return the ``k``-rim of ``ptn`` as a list of integer coordinates.

    The `k`-rim of a partition is the "line between" (or "intersection of") the `k`-boundary and the `k`-interior.  (Section 2.3 of [genocchi]_)

    It will be output as an ordered list of integer coordinates, where the origin is `(0, 0)`.  It will start at the top-left of the `k`-rim (using French convention) and end at the bottom-right.

    EXAMPLES:

    Consider the partition (3, 1) split up into its 1-interior and 1-boundary:

    .. image:: _static/k-rim.JPG
        :height: 180px
        :align: center

    The line shown in bold is the 1-rim, and that information is equivalent to the integer coordinates of the points that occur along that line::

        sage: k_rim(Partition([3, 1]), 1)
        [(3,0), (2,0), (2,1), (1,1), (0,1), (0,2)])

    ..  SEEALSO::

        :meth:`k_interior`, :meth:`k_boundary`, :meth:`boundary`
    """
    interior_rim = boundary(ptn.k_interior(k))
    # get leftmost vertical line
    interior_top_left_y = interior_rim[-1][1]
    v_piece = [(0, y) for y in range(interior_top_left_y + 1, len(ptn) + 1)]
    # get bottommost horizontal line
    interior_bottom_right_x = interior_rim[0][0]
    if ptn:
        ptn_bottom_right_x = ptn[0]
    else:
        ptn_bottom_right_x = 0
    h_piece = [(x, 0)
               for x in range(ptn_bottom_right_x, interior_bottom_right_x, -1)]
    # glue together with boundary
    rim = h_piece + interior_rim + v_piece
    return rim


def k_row_lengths(ptn, k):
    r""" Given a partition, return it's `k`-row-shape.

    This is equivalent to taking the `k`-boundary of the partition and then returning the row-shape of that.  We do *not* discard rows of length 0.  (Section 2.2 of [mem]_)

    EXAMPLES::

        sage: k_row_lengths(Partition([6, 1]), 2)
        [2, 1]

        sage: k_row_lengths(Partition([4, 4, 4, 3, 2]), 2)
        [0, 1, 1, 1, 2]

    ..  SEEALSO::

        :meth:`k_column_lengths`, :meth:`k_boundary`, :meth:`SkewPartition.row_lengths`, :meth:`SkewPartition.column_lengths`
    """
    return ptn.k_boundary(k).row_lengths()


def k_column_lengths(ptn, k):
    r""" Given a partition, return it's `k`-column-shape.

    This is the 'column' analog of :meth:`k_row_lengths`.

    EXAMPLES::

        sage: k_column_lengths(Partition([6, 1]), 2)
        [1, 0, 0, 0, 1, 1]

        sage: k_column_lengths(Partition([4, 4, 4, 3, 2]), 2)
        [1, 1, 1, 2]

    ..  SEEALSO::

        :meth:`k_row_lengths`, :meth:`k_boundary`, :meth:`SkewPartition.row_lengths`, :meth:`SkewPartition.column_lengths`
    """
    return ptn.k_boundary(k).column_lengths()


def has_rectangle(ptn, h, w):
    r""" A partition ``ptn`` has an `h` x `w` rectangle if it's Ferrer's diagram has `h` (*or more*) rows of length `w` (*exactly*).

    EXAMPLES::

        sage: has_rectangle([3, 3, 3, 3], 2, 3)
        True
        sage: has_rectangle([3, 3], 2, 3)
        True
        sage: has_rectangle([4, 3], 2, 3)
        False
        sage: has_rectangle([3], 2, 3)
        False

    ..  SEEALSO::

        :meth:`has_k_rectangle`
    """
    assert h >= 1
    assert w >= 1
    num_rows_of_len_w = 0
    for part in ptn:
        if part == w:
            num_rows_of_len_w += 1
    return num_rows_of_len_w >= h


def has_k_rectangle(ptn, k):
    r""" A partition ``ptn`` has a `k`-rectangle if it's Ferrer's diagram contains `k-i+1` rows (*or more*) of length `i` (*exactly*) for any `i` in `[1, k]`.

    This is mainly a helper function for :meth:`is_k_reducible` and :meth:`is_k_irreducible`, the only difference between this function and :meth:`is_k_reducible` being that this function allows any partition as input while :meth:`is_k_reducible` requires the input to be `k`-bounded.

    EXAMPLES:

    The partition [1, 1, 1] has at least 2 rows of length 1::

        sage: is_k_reducible(Partition([1, 1, 1]), 2)
        True

    The partition [1, 1, 1] does *not* have 4 rows of length 1, 3 rows of length 2, 2 rows of length 3, nor 1 row of length 4::

        sage: is_k_reducible(Partition([1, 1, 1]), 4)
        False

    .. SEEALSO::

        :meth:`is_k_irreducible`, :meth:`is_k_reducible`, :meth:`has_rectangle`
    """
    return any(has_rectangle(ptn, a, b) for (a, b) in k_rectangle_dimension_list(k))


def is_k_bounded(ptn, k):
    r""" Returns ``True`` if and only if the partition ``ptn`` is bounded by ``k``.

    EXAMPLES::

        sage: is_k_bounded(Partition([4, 3, 1]), 4)
        True
        sage: is_k_bounded(Partition([4, 3, 1]), 7)
        True
        sage: is_k_bounded(Partition([4, 3, 1]), 3)
        False
    """
    if ptn.is_empty():
        least_upper_bound = 0
    else:
        least_upper_bound = max(ptn)
    return least_upper_bound <= k


def is_k_reducible(ptn, k):
    r""" A `k`-bounded partition is `k`-*reducible* if it's Ferrer's diagram contains `k-i+1` rows (or more) of length `i` (exactly) for some `i \in [1, k]`.

    (Also, a `k`-bounded partition is `k`-reducible if and only if it is not `k`-irreducible.)

    EXAMPLES:

    The partition [1, 1, 1] has at least 2 rows of length 1::

        sage: is_k_reducible(Partition([1, 1, 1]), 2)
        True

    The partition [1, 1, 1] does *not* have 4 rows of length 1, 3 rows of length 2, 2 rows of length 3, nor 1 row of length 4::

        sage: is_k_reducible(Partition([1, 1, 1]), 4)
        False

    .. SEEALSO::

        :meth:`is_k_irreducible`, :meth:`has_k_rectangle`
    """
    # We only talk about k-reducible / k-irreducible for k-bounded partitions.
    assert is_k_bounded(ptn, k)
    return has_k_rectangle(ptn, k)


def is_k_irreducible(ptn, k):
    r""" A `k`-bounded partition is `k`-*irreducible* if it's Ferrer's diagram does *not* contain `k-i+1` rows (or more) of length `i` (exactly) for every `i \in [1, k]`.

    (Also, a `k`-bounded partition is `k`-irreducible if and only if it is not `k`-reducible.)

    EXAMPLES:

    The partition [1, 1, 1] has at least 2 rows of length 1::

        sage: is_k_irreducible(Partition([1, 1, 1]), 2)
        False

    The partition [1, 1, 1] does *not* have 4 rows of length 1, 3 rows of length 2, 2 rows of length 3, nor 1 row of length 4::

        sage: is_k_irreducible(Partition([1, 1, 1]), 2)
        True

    .. SEEALSO::

        :meth:`is_k_reducible`, :meth:`has_k_rectangle`
    """
    return not is_k_reducible(ptn, k)


def is_symmetric(ptn):
    r"""Given a partition ``ptn``, detect if ``ptn`` equals its own transpose.

    EXAMPLES::

        sage: is_symmetric(Partition([2, 1]))
        True
        sage: is_symmetric(Partition([3, 1]))
        False
    """
    return ptn == ptn.conjugate()


def next_within_bounds(p, min=[], max=None, type=None):
    r"""Get the next partition lexicographically that contains min and is contained in max.

    INPUTS:

    - ``p`` -- The Partition.

    - ``min`` -- (default ``[]``, the empty partition) The 'minimum partition' that ``next_within_bounds(p)`` must contain.

    - ``max`` -- (default ``None``) The 'maximum partition' that ``next_within_bounds(p)`` must be contained in.  If set to ``None``, then there is no restriction.

    - ``type`` -- (default ``None``) The type of partitions allowed.  For example, 'strict' for strictly decreasing partitions, or ``None`` to allow any valid partition.

    EXAMPLES::

        sage: m = [1, 1]
        sage: M = [3, 2, 1]
        sage: Partition([1, 1]).next_within_bounds(min=m, max=M)
        sage: [1, 1, 1]
        sage: Partition([1, 1, 1]).next_within_bounds(min=m, max=M)
        sage: [2, 1]
        sage: Partition([2, 1]).next_within_bounds(min=m, max=M)
        sage: [2, 1, 1]
        sage: Partition([2, 1, 1]).next_within_bounds(min=m, max=M)
        sage: [2, 2]
        sage: Partition([2, 2]).next_within_bounds(min=m, max=M)
        sage: [2, 2, 1]
        sage: Partition([2, 2, 1]).next_within_bounds(min=m, max=M)
        sage: [3, 1]
        sage: Partition([3, 1]).next_within_bounds(min=m, max=M)
        sage: [3, 1, 1]
        sage: Partition([3, 1, 1]).next_within_bounds(min=m, max=M)
        sage: [3, 2]
        sage: Partition([3, 2]).next_within_bounds(min=m, max=M)
        sage: [3, 2, 1]
        sage: Partition([3, 2, 1]).next_within_bounds(min=m, max=M) == None
        sage: True

    ..  SEEALSO::

        :meth:`next`
    """
    # validate inputs
    try:
        assert isinstance(min, (list, Partition))
    except AssertionError:
        raise ValueError('Input parameter ``min`` must be a Partition or a list.')
    try:
        assert isinstance(max, (list, Partition)) or max is None
    except AssertionError:
        raise ValueError('Input parameter ``max`` must be a Partition, a list, or ``None``.')
    # make sure min <= p <= max
    if max is not None:
        assert Partition(max).contains(Partition(p))
    assert Partition(p).contains(Partition(min))
    # check for empty max
    if max is not None and Partition(max).is_empty():
        return None
    # convert partitions to lists to make them mutable
    p = list(p)
    min = list(min)
    # if there is no max, the next partition just tacks a '1' on to the end!
    if max is None:
        return Partition(p + [1])
    # extend p and min to include 0's at the end
    p = p + [0] * (len(max) - len(p))
    min = min + [0] * (len(max) - len(min))
    # finally, run the algo to find next_p
    next_p = copy(p)

    def condition(a, b):
        if type in ('strict', 'strictly decreasing'):
            return a < b - 1
        elif type in (None, 'weak', 'weakly decreasing'):
            return a < b
        else:
            raise ValueError('Unrecognized partition type.')
    for r in range(len(p) - 1, -1, -1):
        if r == 0:
            if (max is None or p[r] < max[r]):
                next_p[r] += 1
                break
            else:
                return None
        else:
            if (max is None or p[r] < max[r]) and condition(p[r], p[r-1]):
                next_p[r] += 1
                break
            else:
                next_p[r] = min[r]
                continue
    return Partition(next_p)


def is_k_core(ptn, k):
    r""" Returns a boolean saying whether or not the Partition ``ptn`` is a ``k``-core.

    EXAMPLES:

    In the partition (2, 1), a hook length of 2 does not occur, but a hook length of 3 does::

        sage: is_k_core(Partition([2, 1]), 2)
        True
        sage: is_k_core(Partition([2, 1]), 3)
        False

    ..  SEEALSO::

        :meth:`Core`
    """
    ptn = Partition(ptn)
    for row_hook_lengths in ptn.hook_lengths():
        for hook_length in row_hook_lengths:
            if hook_length == k:
                return False
    return True


def to_k_core(ptn, k):
    r""" Shift the rows of ``ptn`` minimally in order to create a `k`-core.

    Returns a :class:`Partition` object, not a :class:`Core` object.

    If you plug a `k`-bounded partition into this function and use `k+1` as the input constant, then this is the well-known bijection between `k`-bounded partitions and `k+1`-cores.

    EXAMPLES::

        sage: to_k_core([1, 1], 3)
        [1, 1]
        sage: to_k_core([2, 1], 3)
        [3, 1]
    """
    error = ValueError(
        'The minimal-row-shifting algorithm applied to the partition {} does not produce a {}-core.'.format(ptn, k))
    core = []
    for part in reversed(ptn):
        if core == []:
            core.insert(0, part)
        else:
            core_ptn = Partition(core)
            last_hook_lengths = core_ptn.hook_lengths()[0]
            # this loop could be done away with to make the program more efficient.  You can actually calculate the correct shift amount by looking for k in new_hook_lengths and seeing how much you need to shift.
            minimum_shift = part - previous_part
            for shift in range(minimum_shift, k):
                # the 'shift' is the amount past core[0]
                new_hook_lengths = [l+1+shift for l in last_hook_lengths]
                if k not in new_hook_lengths:
                    # add the appropriate part to core
                    new_part = core[0] + shift
                    # we could improve performance by simply reversing core at the end instead of prepending to the list
                    core.insert(0, new_part)
                    break
            if len(core) == previous_core_len:
                # if none of the shifts were good
                # i think this situation actually can never happen, so if the error occurs, this is a big red flag
                raise error
        previous_part = part
        previous_core_len = len(core)
    core = Partition(core)
    if not is_k_core(core, k):
        raise error
    return core
