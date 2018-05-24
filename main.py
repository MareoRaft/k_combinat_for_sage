#!/usr/bin/env sage
from sage.all import *

# HELPERS (only exist as helper functions for other things):
def is_weakly_decreasing(li):
    return all(li[i] >= li[i+1] for i in range(len(li)-1))


# MAIN:


# SkewPartition methods:
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


# Partition methods:
def k_row_lengths(ptn, k):
    return ptn.k_boundary(k).row_lengths()

def k_column_lengths(ptn, k):
    return ptn.k_boundary(k).column_lengths()

def has_rectangle(ptn, h, w):
    """ A partition has an `h` x `w` rectangle if it's Ferrer's diagram has `h` (*or more*) rows of length `w` (*exactly*).
    """
    assert h >= 1
    assert w >= 1
    num_rows_of_len_w = 0
    for part in ptn:
        if part == w:
            num_rows_of_len_w += 1
    return num_rows_of_len_w >= h

def has_k_rectangle(ptn, k):
    """ A partition has a k-rectangle if it's Ferrer's diagram contains k-i+1 rows (or more) of length i (exactly) for any i in [1, k].
    """
    for i in range(1, k+1):
        if has_rectangle(ptn, k-i+1, i):
            return True
    return False

def is_k_bounded(ptn, k):
    """ Returns True iff the partition is bounded by k. """
    if ptn.is_empty():
        least_upper_bound = 0
    else:
        least_upper_bound = max(ptn)
    return least_upper_bound <= k

def Partition_is_k_reducible(ptn, k):
    """ A k-bounded partition is __k-reducible__ if it has a k-rectangle. """
    # We only talk about k-reducible / k-irreducible for k-bounded partitions.
    assert is_k_bounded(ptn, k)
    return has_k_rectangle(ptn, k)

def Partition_is_k_irreducible(ptn, k):
    return not Partition_is_k_reducible(ptn, k)
# END Partition methods.

def get_k_rectangles(k):
    """ A __k-rectangle__ is a partition whose Ferrer's diagram is a rectangle whose largest hook-length is k. """
    return [Partition([i] * (k-i+1)) for i in range(1, k+1)]

def get_k_irreducible_partition_lists(k):
    """Since there are n! such partitions, the big-O time can't be better than that.
    We could have a yeild in the function to be an iterator.

    Returns: list of lists (instead of list of Partition objects)
    """
    k = NN(k)
    k_irr_ptns = [[]]
    # NO rows of length k
    for i in range(1, k):
        new_k_irr_ptns = []
        for ptn in k_irr_ptns:
            # at most i rows of length k-i where 1 <= i < k
            for num_rows in range(0, i+1):
                new_ptn = ptn + [k-i]*num_rows
                new_k_irr_ptns.append(new_ptn)
        k_irr_ptns = new_k_irr_ptns
    return k_irr_ptns

def get_k_irreducible_partitions(k):
    """Given k, return the n! k-irreducible-partitions. """
    return [Partition(e) for e in get_k_irreducible_partition_lists(k)]

def n_to_number_of_linked_partition_self_pairs(n):
    """ Given a natural number n, count how many partitions l of size n have the property that (l, l) has a corresponding linked-skew-diagram. """
    ps = Partitions(n)
    count = 0
    for p in ps:
        try:
            row_col_to_skew_partition(p, p)
        except:
            pass
        else:
            count += 1
    return count

def bump_path_piece(sp, start_row_index, blocked_rows=set()):
    # this algo find the correct "L" piece of the path, where the bottom right cell is cell1, the bottom left is cell2, and the top left is cell3
    # Returns (top_row_index, is_end) which are the row index of cell3 and whether or not we 'broke free' out of the top left cell of the skew-partition, respectively.
    col_index2 = left(sp, start_row_index)
    row_index3 = top(sp, col_index2) + 1
    while row_index3 in blocked_rows:
        row_index3 += 1
    # CATTY-CORNER ONLY line:
    max_row_index = len(sp.outer()) - 1
    if row_index3 > max_row_index:
        return None, True
    else:
        return row_index3, False
def bump_path(sp, row_index, blocked_rows=set()):
    new_blocked_rows = {row_index}
    # is_end says if you reached the end of the path
    while True:
        row_index, is_end = bump_path_piece(sp, row_index, blocked_rows)
        if is_end:
            break
        else:
            new_blocked_rows.add(row_index)
    return new_blocked_rows
def skew_partition_to_selected_rows(sp):
    # actually this may ONLY WORK for catty-connected skew-partitions, because i'm not sure how we deal with 'missing' rows
    # arguably we should call it a linked_skew_partition
    # record the indices of rows that have been used up
    selected_rows = set()
    blocked_rows = set()
    for row_index, outer_row in enumerate(sp.outer()):
        if row_index not in blocked_rows:
            selected_rows.add(row_index)
            new_blocked_rows = bump_path(sp, row_index, blocked_rows)
            blocked_rows.update(new_blocked_rows)
    return sorted(selected_rows)
def selected_rows_to_root_ideal(n, selected_indecis):
    """Given the dimension of the square n and the selected rows, output the root ideal """
    root_ideal_cells = []
    selected_indecis = set(selected_indecis)
    permitted_col_indecis = set(range(n)) - selected_indecis
    for i in range(n):
        if i in selected_indecis:
            if permitted_col_indecis:
                smallest_unblocked_index = min(permitted_col_indecis)
                root_ideal_cells += [(i, j) for j in range(smallest_unblocked_index, n)]
                permitted_col_indecis.remove(smallest_unblocked_index)
                selected_indecis.add(smallest_unblocked_index)
    return root_ideal_cells
def skew_partition_to_root_ideal(sp):
    selected_indecis = skew_partition_to_selected_rows(sp)
    n = len(sp.outer())
    root_ideal = selected_rows_to_root_ideal(n, selected_indecis)
    return root_ideal

def is_symmetric(l):
    """Given a partition l, detect if l = l'.

    This function runs in LINEAR time of order length(l).
    """
    for j in range(0, len(l)):
        for k in range(l[-j], l[-j-1]):
            if l[k] != len(l) - j:
                return False
    return True

def sequence(func, num_terms=20):
    seq = []
    for n in range(0, num_terms):
        seq.append(func(n))
    return seq

def print_sequence(func, num_terms=float('inf')):
    n = 0
    while n < num_terms:
        print('n={}\t{}=f(n)'.format(n, func(n)))

def is_k_shape(ptn, k):
    """ A partition is a k-shape if its k-boundary has row-shape and col-shape that are partitions themselves. """
    if k is None:
        # see if it's a k-shape for any k in [1, n-1].
        # (note that every partition is a 0-shape and an n-shape)
        n = ptn.size()
        lis = [is_k_shape(ptn, kk) for kk in range(1, n)]
        return any(lis)
    else:
        k_bdy = ptn.k_boundary(k)
        return is_linked(k_bdy)

def n_to_k_shapes(n, k):
    """ Given n, find all partitions of size n that are k-shapes. """
    return [ptn for ptn in Partitions(n) if is_k_shape(ptn, k)]

def n_to_num_k_shapes(n, k):
    return len(n_to_k_shapes(n, k))

def n_to_k_skews(n, k):
    """ Given n, find all k-skews coming from partitions of size n. """
    return [ptn.k_boundary(k) for ptn in Partitions(n) if is_k_shape(ptn, k)]

def n_to_self_conjugate_k_skews(n, k):
    k_skews = n_to_k_skews(n, k)
    return [ks for ks in k_skews if ks == ks.conjugate()]

def n_to_num_self_conjugate_k_skews(n, k):
    return len(n_to_self_conjugate_k_skews(n, k))



# kShape methods:
def kShape_is_k_reducible(s, k):
    """ A k-shape is called __k-reducible__ if there exists a k-rectangle R such that (the k-row-shape has R and the k-column-shape has R'). """
    def has_k_rectangle_pair(k):
        for i in range(1, k+1):
            a = k-i+1
            b = i
            if has_rectangle(rs, a, b) and has_rectangle(cs, b, a):
                return True
        return False
    rs = Partition(k_row_lengths(s, k))
    cs = Partition(k_column_lengths(s, k))
    return has_k_rectangle_pair(k) or has_k_rectangle_pair(k-1)

def kShape_is_k_irreducible(s, k):
    """ A k-shape is called __k-irreducible__ if it is not k-reducible. """
    return not kShape_is_k_reducible(s, k)
# END k-shape methods.

def get_k_irreducible_k_shapes(k):
    # The k-row-shape has at most k rows of length 0, k-1 rows of length 1, ..., 0 rows of length k.  And 0 rows of length greater than k.  Hence the k-row-shape has an upper bound of k*(k-1)/2 rows.  The same goes for the k-col-shape.
    bound = (k-1)*k/2
    n_bound = bound**2
    ptns = []
    for n in range(0, n_bound+1):
        ptns += Partitions(n, max_length=bound, max_part=bound)
        k_irr_k_shapes = [p for p in ptns if is_k_shape(p, k) and kShape_is_k_irreducible(p, k)]
    return k_irr_k_shapes

