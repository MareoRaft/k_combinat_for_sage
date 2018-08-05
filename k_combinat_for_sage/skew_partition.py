# -*- coding: utf-8 -*-
r"""
Sage has a builtin `SkewPartition <https://doc.sagemath.org/html/en/reference/combinat/sage/combinat/skew_partition.html>`_ object (a.k.a. skew-shape).  *This* module adds extra useful functions for skew partitions:

REFERENCES:

.. [mem] Lam, T., Lapointe, L., Morse, J., & Shimozono, M. (2013). `The poset of k-shapes and branching rules for k-Schur functions <http://breakfreerun.org/index.php/ebooks/the-poset-of-k-shapes-and-branching-rules-for-k-schur-functions>`_. Memoirs of the American Mathematical Society, 223(1050), 1-113. DOI: 10.1090/S0065-9266-2012-00655-1
"""
from sage.all import *
# ^*^ sphinx insert ^*^

# HELPERS (only exist as helper functions for other things):


def is_weakly_decreasing(li):
    return all(li[i] >= li[i+1] for i in range(len(li)-1))


def is_strictly_decreasing(li):
    return all(li[i] > li[i+1] for i in range(len(li)-1))


# SkewPartition stuff
def is_symmetric(sp):
    r""" A SkewPartition is *symmetric* if its inner and outer shapes are symmetric.

    Returns True if and only if the SkewPartition `sp` is equal to its own conjugate.
    """
    return sp == sp.conjugate()


def right(sp, row_index):
    r""" Given a SkewPartition and a 0-based row index, return the 0-based column index of the *rightmost* cell in the corresponding row.  (Section 2.1 of [mem]_)

    EXAMPLES::

        sage: right(SkewPartition([[4, 1], [2]]), 0)
        3

        sage: right(SkewPartition([[4, 1], [2]]), 1)
        0

    An input index that is out of the bounds of the skew partition will throw an error::

        sage: right(SkewPartition([[4, 1], [2]]), 2)
        IndexError: list index out of range

    An in-bounds index where no cells exist will return ``None``::

        sage: right(SkewPartition([[2, 1, 1, 1], [1, 1]]), 1) == None
        True
    """
    # first check to make sure the cell exists
    if sp.row_lengths()[row_index] == 0:
        return None
    outer_row_lengths = sp.outer().to_list()
    outer_row_length = outer_row_lengths[row_index]
    col_index = outer_row_length - 1
    return col_index


def left(sp, row_index):
    r""" Given a SkewPartition and a 0-based row index, return the 0-based column index of the *leftmost* cell in the corresponding row.  (Section 2.1 of [mem]_)

    EXAMPLES::

        sage: left(SkewPartition([[4, 1], [2]]), 0)
        2

        sage: left(SkewPartition([[4, 1], [2]]), 1)
        0

    An input index that is out of the bounds of the skew partition will throw an error::

        sage: left(SkewPartition([[4, 1], [2]]), 2)
        IndexError: list index out of range

    An in-bounds index where no cells exist will return ``None``::

        sage: left(SkewPartition([[2, 1, 1, 1], [1, 1]]), 1) == None
        True
    """
    # first check to make sure the cell exists
    if sp.row_lengths()[row_index] == 0:
        return None
    outer_row_lengths = sp.outer().to_list()
    inner_row_lengths = sp.inner().to_list()
    if (row_index in range(len(outer_row_lengths))) and (row_index not in range(len(inner_row_lengths))):
        inner_row_length = 0
    else:
        inner_row_length = inner_row_lengths[row_index]
    col_index = inner_row_length
    return col_index


def top(sp, col_index):
    r""" Given a SkewPartition and a 0-based column index, return the 0-based row index of the *topmost* cell in the corresponding column.  (Section 2.1 of [mem]_)

    EXAMPLES::

        sage: top(SkewPartition([[4, 1, 1], [2]]), 0)
        2

    An in-bounds col_index where no cells exist will return ``None``::

        sage: top(SkewPartition([[4, 1, 1], [2]]), 1) == None
        True

        sage: top(SkewPartition([[4, 1, 1], [2]]), 2)
        0

        sage: top(SkewPartition([[4, 1, 1], [2]]), 3)
        0

    A col_index that is out-of-bounds of the skew partition will throw an error::

        sage: top(SkewPartition([[4, 1, 1], [2]]), 4)
        IndexError: list index out of range
    """
    return right(sp.conjugate(), col_index)


def bottom(sp, col_index):
    r""" Given a SkewPartition and a 0-based column index, return the 0-based row index of the *bottommost* cell in the corresponding column.  (Section 2.1 of [mem]_)

    EXAMPLES::

        sage: bottom(SkewPartition([[4, 1, 1], [2]]), 0)
        1

    An in-bounds col_index where no cells exist will return ``None``::

        sage: bottom(SkewPartition([[4, 1, 1], [2]]), 1) == None
        True

        sage: bottom(SkewPartition([[4, 1, 1], [2]]), 2)
        0

        sage: bottom(SkewPartition([[4, 1, 1], [2]]), 3)
        0

    A col_index that is out-of-bounds of the skew partition will throw an error::

        sage: bottom(SkewPartition([[4, 1, 1], [2]]), 4)
        IndexError: list index out of range
    """
    return left(sp.conjugate(), col_index)


def is_linked(sp):
    r"""
    A skew-shape ``sp`` is a *skew-linked diagram* if both the row-shape and column-shape of `sp` are partitions.

    EXAMPLES:

    Both row shape and column shape are valid::

        sage: is_linked(SkewPartition([[2, 1], [1]]))
        True

    Valid row shape but invalid column shape::

        sage: is_linked(SkewPartition([[3, 2], [1]]))
        False
    """
    return is_weakly_decreasing(sp.row_lengths()) and is_weakly_decreasing(sp.column_lengths())


def row_col_to_skew_partition(rs, cs):
    # this ALREADY exists in sage. see SkewPartition.from_row_and_column_length
    outer = []
    inner = []
    current_cs = [0] * len(cs)
    row_index = 0
    for col_coindex, col_length in enumerate(list(reversed(cs))):
        current_col_length = list(reversed(current_cs))[col_coindex]
        num_rows_to_slide = col_length - current_col_length
        if num_rows_to_slide < 0:
            raise ValueError(
                'The inputted (row-shape, col-shape) pair has no possible corresponding skew-shape.')
        # 'col_num' is 1-based index of cols
        col_num = len(cs) - col_coindex
        while num_rows_to_slide > 0:
            if row_index > len(rs) - 1:
                raise ValueError('error more')
            # slide a row
            outer.append(col_num)
            inner.append(col_num - rs[row_index])
            # update params/info
            for c in range(col_num - rs[row_index], col_num):
                current_cs[c] += 1
            row_index += 1
            num_rows_to_slide -= 1
    return SkewPartition([outer, inner])


def k_boundary_to_partition(sp, k=None, strict=True):
    r""" Given a ``k``-boundary ``sp`` (`k`-boundaries are a specific type of skew-shape), output the original partition whose `k`-boundary is `sp`.

    (For the definition of `k`-boundary, see Section 2.2 of [mem]_)

    If strict is set to True, the program will assert that the skew-shape really is a `k`-boundary.

    TODO: test

    EXAMPLES::

        sage: k_boundary_to_partition(SkewPartition([[3, 2, 1], [2, 1]]))
        [3, 2, 1]

        sage: k_boundary_to_partition(SkewPartition([[3, 1], [2]]), k=2)
        Error
        sage: k_boundary_to_partition(SkewPartition([[3, 1], [2]]), strict=False)
        [3, 1]
    """
    if strict:
        assert is_k_boundary(sp, k)
    return sp.outer()


def is_k_boundary(sp, k=None):
    r""" Given a skew-shape ``sp`` and natural number ``k``, return True if and only if `sp` is a `k`-boundary.  (Section 2.2 of [mem]_)

    Given a skew-shape `sp` *only*, return True if and only if there exists some `k` such that `sp` is a `k`-boundary.

    TODO: test

    EXAMPLES::

        sage: is_k_boundary(SkewPartition([[3, 2, 1], [2, 1]]))
        True
        sage: is_k_boundary(SkewPartition([[3, 2, 1], [2, 1]]), k=1)
        True
        sage: is_k_boundary(SkewPartition([[3, 2, 1], [2, 1]]), k=2)
        False
    """
    if k is None:
        max_hook_length = sp.outer().hook_length(0, 0)
        return any(is_k_boundary(sp, k_star) for k_star in range(0, max_hook_length+1))
    elif k == 0:
        # the only valid 0-boundary is the empty shape
        return sp.outer() == sp.inner()
    else:
        r"""We go down and left of each cell to create the only possible partition that could have led to this potential k-boundary

        (Any other partition containing this skew_shape would necessarily have a northeast corner that the skew_shape does *not have*.  But in order for the skew-shape to be a k-boundary, it *must have* that northeast corner.)
        """
        l = k_boundary_to_partition(sp, strict=False)
        r"""now that we have the partition, we simply compute it's hook-length for each cell and verify that for each cell of values k or less, it appears in the sp"""
        correct_k_boundary = l.k_boundary(k)
        return sp == correct_k_boundary


def row_shape_to_linked_skew_partitions(rs):
    r""" Given a partition ``rs``, find all linked SkewPartitions whose row-shape is ``rs``.

    EXAMPLES:

    Note that [4, 2, 1] / [1, 1] is *not* linked and hence doesn't appear in the list below::

        sage: row_shape_to_linked_skew_partitions(Partition([3, 1, 1]))
        [[3, 1, 1] / [], [4, 1, 1] / [1], [5, 2, 1] / [2, 1]]
    """
    def ptn_to_linked_things(p):
        def thing_to_added_row_things(sp, row_len):
            def add_row(sp, row_len, offset):
                # add the next row onto the skew shape
                outer = sp.outer().to_list()
                inner = sp.inner().to_list()
                inner += [0] * (len(outer) - len(inner))
                outer = [e + offset for e in outer]
                inner = [e + offset for e in inner]
                outer.append(row_len)
                return SkewPartition([outer, inner])
            # START thing_to_added_row_things
            previous_checked_col_index = sp.outer()[-1]
            # find the maximum leftmost offset for the new row
            col_lens = sp.column_lengths()
            max_offset = row_len
            prev_col_len = 0
            for col_index in range(previous_checked_col_index, -1, -1):
                # get length of column
                col_len = col_lens[col_index] if col_index < len(col_lens) else 0
                # check col-shape partition condition
                if col_len >= prev_col_len:
                    # col_index is good, continue
                    prev_col_len = col_len
                else:
                    # col_index is bad, stop
                    good_col_index = col_index + 1
                    max_offset = row_len - good_col_index
                    break
            # now add all possible positions for the row onto the list
            return [add_row(sp, row_len, offset) for offset in range(0, max_offset+1)]
        # START ptn_to_linked_things
        assert isinstance(p, list) and not isinstance(p, Partition)
        if len(p) <= 1:
            return [SkewPartition([p, []])]
        else:
            # these incomplete guys are not necessarily skew partitions
            incomplete_things = ptn_to_linked_things(p[:-1])
            almost_complete_things = []
            for incomplete_thing in incomplete_things:
                almost_complete_things += thing_to_added_row_things(
                    incomplete_thing, p[-1])
            return almost_complete_things
    # START row_shape_to_linked_skew_partitions
    rs_zero = list(Partition(rs)) + [0]
    return ptn_to_linked_things(rs_zero)


def size_to_linked_skew_partitions(size):
    r""" Given a natural number ``size``, return all linked SkewPartitions of size ``size``.

    EXAMPLES::

        sage: size_to_linked_skew_partitions(3)
        [[3] / [], [2, 1] / [], [3, 1] / [1], [1, 1, 1] / [], [2, 1, 1] / [1], [3, 2, 1] / [2, 1]]
    """
    linked_skew_ptns = []
    # Here is one major optimization that's possible: Instead of first calculating all Partitions(size), and then doing the ptn_to_linked_things algo for each partition, actually go through the work of generating the partitions manually, and use ptn_to_linked_things algo as you go.  This is to ELIMINATE the redundancy of having two partitions that START with the same sub-partition.
    ptns = Partitions(size)
    for ptn in ptns:
        linked_skew_ptns += row_shape_to_linked_skew_partitions(ptn)
    return linked_skew_ptns
