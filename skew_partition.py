# -*- coding: utf-8 -*-
from sage.all import *

# HELPERS (only exist as helper functions for other things):
def is_weakly_decreasing(li):
    return all(li[i] >= li[i+1] for i in range(len(li)-1))

def is_strictly_decreasing(li):
    return all(li[i] > li[i+1] for i in range(len(li)-1))


# SkewPartition stuff
def right(self, row_index):
    """self: a SkewPartition
    Given a 0-based row_index, return the 0-based column index of the rightmost cell in the corresponding row """
    # first check to make sure the cell exists
    if self.row_lengths()[row_index] == 0:
        return None
    outer_row_lengths = self.outer().to_list()
    outer_row_length = outer_row_lengths[row_index]
    col_index = outer_row_length - 1
    return col_index

def left(self, row_index):
    """self: a SkewPartition
    Given a 0-based row_index, return the 0-based column index of the leftmost cell in the corresponding row """
    # first check to make sure the cell exists
    if self.row_lengths()[row_index] == 0:
        return None
    outer_row_lengths = self.outer().to_list()
    inner_row_lengths = self.inner().to_list()
    if (row_index in range(len(outer_row_lengths))) and (row_index not in range(len(inner_row_lengths))):
        inner_row_length = 0
    else:
        inner_row_length = inner_row_lengths[row_index]
    col_index = inner_row_length
    return col_index

def top(self, col_index):
    """self: a SkewPartition
    Given a 0-based col_index, return the 0-based row_index of the topmost cell in the corresponding column """
    return right(self.conjugate(), col_index)

def bottom(self, col_index):
    """self: a SkewPartition
    Given a 0-based col_index, return the 0-based row_index of the bottommost cell in the corresponding column """
    return left(self.conjugate(), col_index)

def is_linked(self):
    """
    A skew-shape `s` is a __skew-linked diagram__ if both the row-shape and column-shape of `s` are partitions.
    """
    return is_weakly_decreasing(self.row_lengths()) and is_weakly_decreasing(self.column_lengths())

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
            raise ValueError('The inputted (row-shape, col-shape) pair has no possible corresponding skew-shape.')
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



"""
Given a skew-linked diagram, is it a k-boundary?  (That is, does there exist some partition which - when cells of hook-length > k are removed - becomes the skew-linked diagram.)
"""

def k_boundary_to_partition(skew_shape, strict=True):
    """
    skew_shape: The skew-shape (a k-boundary) to find the corresponding partition for.
    strict: If true, assert that the skew-shape really is a k-boundary.
    """
    if strict:
        assert is_k_boundary(skew_shape)
    # return minimal_containing_partition()
    return skew_shape.outer()

def is_k_boundary(skew_shape, k):
    if k == 0:
        # the only valid 0-boundary is the empty shape
        return skew_shape.outer() == skew_shape.inner()
    else:
        """We go down and left of each cell to create the only possible partition that could have led to this potential k-boundary

        (Any other partition containing this skew_shape would necessarily have a northeast corner that the skew_shape does *not have*.  But in order for the skew-shape to be a k-boundary, it *must have* that northeast corner.)
        """
        l = k_boundary_to_partition(skew_shape, strict=False)
        """now that we have the partition, we simply compute it's hook-length for each cell and verify that for each cell of values k or less, it appears in the skew_shape"""
        correct_k_boundary = l.k_boundary(k)
        return skew_shape == correct_k_boundary





########### GETTER FUNCS ##############
def add_row(sp, row_len, offset):
    outer = sp.outer().to_list()
    inner = sp.inner().to_list()
    inner += [0] * (len(outer) - len(inner))
    outer = [e + offset for e in outer]
    inner = [e + offset for e in inner]
    outer.append(row_len)
    return SkewPartition([outer, inner])
def thing_to_added_row_things(sp, row_len):
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
def ptn_to_linked_things(p):
    assert isinstance(p, list) and not isinstance(p, Partition)
    if len(p) <= 1:
        return [SkewPartition([p, []])]
    else:
        # these incomplete guys are not necessarily skew partitions
        incomplete_things = ptn_to_linked_things(p[:-1])
        almost_complete_things = []
        for incomplete_thing in incomplete_things:
            almost_complete_things += thing_to_added_row_things(incomplete_thing, p[-1])
        return almost_complete_things
def ptn_to_linked_skew_partitions(p):
    """ Given a partition p, find all linked SkewPartitions whose row-shape is the partition. """
    p_zero = list(Partition(p)) + [0]
    return ptn_to_linked_things(p_zero)

def n_to_linked_skew_partitions(n):
    """ Given n, return all linked SkewPartitions of size n. """
    linked_skew_ptns = []
    # Here is one major optimization that's possible: Instead of first calculating all Partitions(n), and then doing the ptn_to_linked_things algo for each partition, actually go through the work of generating the partitions manually, and use ptn_to_linked_things algo as you go.  This is to ELIMINATE the redundancy of having two partitions that START with the same sub-partition.
    ptns = Partitions(n)
    for ptn in ptns:
        linked_skew_ptns += ptn_to_linked_skew_partitions(ptn)
    return linked_skew_ptns
