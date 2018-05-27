#!/usr/bin/env sage
# -*- coding: utf-8 -*-
from sage.all import *

# HELPERS (only exist as helper functions for other things):
def is_weakly_decreasing(li):
    return all(li[i] >= li[i+1] for i in range(len(li)-1))

def is_strictly_decreasing(li):
    return all(li[i] > li[i+1] for i in range(len(li)-1))

def k_rectangle_dimension_list(k):
    return [(k-i+1, i) for i in range(1, k+1)]

def get_n_from_root_ideal(root_ideal):
    return max(c for (r,c) in root_ideal) + 1


# MAIN:


class RootIdeal(list):
    def __hash__(self):
        return hash(tuple(sorted(self)))
RI = RootIdeal


# SkewPartition methods:
def SkewPartition_right(self, row_index):
    """self: a SkewPartition
    Given a 0-based row_index, return the 0-based column index of the rightmost cell in the corresponding row """
    # first check to make sure the cell exists
    if self.row_lengths()[row_index] == 0:
        return None
    outer_row_lengths = self.outer().to_list()
    outer_row_length = outer_row_lengths[row_index]
    col_index = outer_row_length - 1
    return col_index

def SkewPartition_left(self, row_index):
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

def SkewPartition_top(self, col_index):
    """self: a SkewPartition
    Given a 0-based col_index, return the 0-based row_index of the topmost cell in the corresponding column """
    return SkewPartition_right(self.conjugate(), col_index)

def SkewPartition_bottom(self, col_index):
    """self: a SkewPartition
    Given a 0-based col_index, return the 0-based row_index of the bottommost cell in the corresponding column """
    return SkewPartition_left(self.conjugate(), col_index)

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

def bump_path_piece(sp, start_row_index, blocked_rows=set()):
    # this algo find the correct "L" piece of the path, where the bottom right cell is cell1, the bottom left is cell2, and the top left is cell3
    # Returns (top_row_index, is_end) which are the row index of cell3 and whether or not we 'broke free' out of the top left cell of the skew-partition, respectively.
    col_index2 = SkewPartition_left(sp, start_row_index)
    row_index3 = SkewPartition_top(sp, col_index2) + 1
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
def selected_rows_to_maximum_root_ideal(n, selected_indecis):
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

def skew_partition_to_removable_roots(sp, type='max'):
    """
    sp: The SkewPartition.
    type: Get removable roots for the 'min' root ideal or the 'max' root ideal.
    returns: A list of removable roots.

    reference: [LM04] L. Lapointe and J. Morse. Order ideals in weak subposets of Young’s lattice and associated
unimodality conjectures. Ann. Comb., 8(2):197–219, 2004.
    """
    def type_shift(x, type):
        if type == 'max':
            return x
        elif type == 'min':
            return x - 1
        else:
            raise ValueError('Bad type.')
    assert is_linked(sp)
    mu = Partition(sp.row_lengths())
    eta = sp.inner()
    return [(i, mu[type_shift(eta[i], type)] + i) for i in range(0, len(eta))]

def removable_roots_to_partition(corners, n):
    corners = sorted(corners)
    # r is the row index or the 'y' value
    # c is the col index of the 'x' value
    previous_r = -1
    ptn = []
    for r in range(0, n):
        if corners:
            current_r = corners[0][0]
            current_c = corners[0][1]
            if current_r == r:
                # see how many rows to fill
                num_rows = r - previous_r
                ptn += [n - current_c] * num_rows
                # delete corner since it has now been used
                previous_r = current_r
                corners = corners[1:]
    return Partition(ptn)

def removable_roots_to_root_ideal(corners, n):
    ptn = removable_roots_to_partition(corners, n)
    ri = partition_to_root_ideal(ptn, n)
    return ri

def skew_partition_to_root_ideal(sp, type='max', method='removable roots'):
    if method == 'removable roots':
        corners = skew_partition_to_removable_roots(sp, type)
        n = len(sp.outer())
        root_ideal = removable_roots_to_root_ideal(corners, n)
    elif method == 'bounce':
        if type != 'max':
            raise Exception('The bounce method can only yield the maximum root ideal (type=max).')
        selected_indecis = skew_partition_to_selected_rows(sp)
        n = len(sp.outer())
        root_ideal = selected_rows_to_maximum_root_ideal(n, selected_indecis)
    else:
        raise ValueError('Unknown method.')
    return RI(root_ideal)

def RootIdeal_next(ri, min=[], max=None, n=None):
    # figure out dimension of square
    if n is not None:
        pass
    elif ri:
        n = get_n_from_root_ideal(ri)
    elif min:
        n = get_n_from_root_ideal(min)
    elif max:
        n = get_n_from_root_ideal(max)
    else:
        raise Exception('There is no way to figure out the size of the staircase that the root ideals fall in.  Please supply n.')
    ptn = root_ideal_to_partition(ri)
    min_ptn = root_ideal_to_partition(min)
    max_ptn = root_ideal_to_partition(max)
    next_ptn = Partition_next(ptn, min=min_ptn, max=max_ptn)
    next_ri = partition_to_root_ideal(next_ptn, n)
    return next_ri

def skew_partition_to_root_ideals(sp):
    """ Given a skew partition, find the corresponding set (but given as a list here) of root_ideals.

    We could change this to an iterator if users may not want all the root ideals.
    """
    min_ri = skew_partition_to_root_ideal(sp, type='min')
    max_ri = skew_partition_to_root_ideal(sp, type='max')
    n = len(sp.outer())
    next_func = lambda ri: RootIdeal_next(ri, min=min_ri, max=max_ri, n=n)
    ris = generate_path(next_func, min_ri)
    return ris

def down(root_ideal, row_index):
    """ Given a root ideal and a starting position (row_index), move right unti you hit the root ideal, then move straight down until you hit the diagonal, and return the new index. """

    # Note: I am assuming the cells in the root ideal are IN ORDER with y coordinates weakly increasing, and for fixed y, x strictly increasing
    for (r,c) in root_ideal:
        if r == row_index:
            return c
    return None

def up(root_ideal, index):
    for (r,c) in reversed(root_ideal):
        if c == index:
            return r
    return None

def generate_path(next_func, start):
    path = [start]
    while True:
        next_ = next_func(path[-1])
        if next_ is not None:
            path.append(next_)
        else:
            break
    return path

def down_path(root_ideal, start_index):
    next_func = lambda index: down(root_ideal, index)
    return generate_path(next_func, start_index)

def up_path(root_ideal, start_index):
    next_func = lambda index: up(root_ideal, index)
    return generate_path(next_func, start_index)

def top(root_ideal, start_index):
    return up_path(root_ideal, start_index)[-1]

def bottom(root_ideal, start_index):
    return down_path(root_ideal, start_index)[-1]

def down_path_partition_sum(root_ideal, ptn, start_index):
    """ This is \\mu_i in Definition 2.3 of SKEW-LINKED CATALAN FUNCTIONS AND k-SCHUR POSITIVITY. """
    return sum(ptn[j] for j in down_path(root_ideal, start_index))
def down_path_partition(root_ideal, ptn):
    """ This is the *column shape* \\mu' as defined by Definition 2.3 of SKEW-LINKED CATALAN FUNCTIONS AND k-SCHUR POSITIVITY.  It is also introduced in the second paragraph of the overview as \\mathfrak{cs}(\Psi, \lambda). """
    if not root_ideal:
        mu = ptn
    else:
        mu = []
        # n is the side length of the square
        n = get_n_from_root_ideal(root_ideal)
        indecis_available = set(range(0, n))
        for index in range(0, n):
            if index in indecis_available:
                # add the kthing to mu
                mu.append(down_path_partition_sum(root_ideal, ptn, index))
                # remove indecis from future draws
                dpath = down_path(root_ideal, index)
                indecis_available -= set(dpath)
    return Partition(mu)

def root_ideal_to_partition(root_ideal):
    """ Given a root ideal (list of cells), return the corresponding partition (the row shape of the root ideal). """
    if root_ideal is None or root_ideal == False:
        return root_ideal
    if not root_ideal:
        ptn = []
    else:
        max_r = root_ideal[-1][0]
        ptn = [0] * (max_r + 1)
        for (r,c) in root_ideal:
            ptn[r] += 1
    return Partition(ptn)

def partition_to_root_ideal(ptn, n):
    """ Given a partition and the size of the square, return the corresponding root ideal.  (This is the inverse function to root_ideal_to_partition when restricted to n x n grid and sub-(n-1)-staircase partitions.) """
    if ptn is None or ptn == False:
        return ptn
    root_ideal = []
    for r, part in enumerate(ptn):
        root_ideal += [(r, c) for c in range(n-part, n)]
    return RI(root_ideal)

def is_rational_root_ideal(ri):
    """ Given a root ideal ri, check to see if it is a *rational root ideal*, as defined in Example 2.4 of SKEW-LINKED CATALAN FUNCTIONS AND k-SCHUR POSITIVITY.  This merely means that it's corresponding partition is strictly decreasing! """
    ptn = root_ideal_to_partition(ri)
    return is_strictly_decreasing(ptn)



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
# End SkewPartition methods.


# Partition methods:
def boundary(ptn):
    """ The boundary of a partition is the set { NE(d) | for all d diagonal }.  That is, for every diagonal, we find the northeasternmost (NE) point on that diagonal which is also in the ferrer's diagram.  Finally, we restrict to only integer coordinates.

    The boundary will go from bottom-right to top-left in the French notation.
    """
    def horizontal_piece((start_x, start_y), bdy):
        if not bdy:
            h_piece = [(start_x, start_y)]
        else:
            stop_x = bdy[-1][0]
            y = start_y # y never changes
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
    # claim victory
    return bdy

def k_rim(ptn, k):
    """ The k-rim of a partition is the "line between" (or "intersection of") the k-boundary and the k-interior.

    It will be outputted as an ordered list of integer coordinates, where the origin is (0,0).  It will start at the top-left of the k-rim (under French depiction) and end at the bottom-right. """
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
    h_piece = [(x, 0) for x in range(ptn_bottom_right_x, interior_bottom_right_x, -1)]
    # glue together with boundary
    rim = h_piece + interior_rim + v_piece
    return rim

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
    return any(has_rectangle(ptn, a, b) for (a, b) in k_rectangle_dimension_list(k))

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

def is_symmetric(ptn):
    """Given a partition λ, detect if λ = λ'.

    This function runs in LINEAR time of order length(λ).
    """
    for j in range(0, len(ptn)):
        for k in range(ptn[-j], ptn[-j-1]):
            if ptn[k] != len(ptn) - j:
                return False
    return True

def Partition_next(p, min=[], max=None):
    """
    Get the next partition lexigraphically that contains min and is contained in max.
    ptn: The Partition.
    min: The 'minimum partition' that next_advanced(ptn) must contain.
    max: The 'maximum partition' that next_advanced(ptn) must be contained in.
    """
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
    for r in range(len(p) - 1, -1, -1):
        if r == 0:
            if (max is None or p[r] < max[r]):
                next_p[r] += 1
                break
            else:
                return None
        else:
            if (max is None or p[r] < max[r]) and p[r] < p[r-1]:
                next_p[r] += 1
                break
            else:
                next_p[r] = min[r]
                continue
    return Partition(next_p)
# END Partition methods.

def get_k_rectangles(k):
    """ A __k-rectangle__ is a partition whose Ferrer's diagram is a rectangle whose largest hook-length is k. """
    return [Partition([a] * b) for (a, b) in k_rectangle_dimension_list(k)]

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
def h_bounds(p, k, width):
    """ Recall the "H_i" as defined in Def 3.3 of Combinatorics of k-shapes and Genocchi
numbers.
    k: The k used for the 'k'-shape or 'k'-boundary
    width: i
    returns: (y_min, y_max) The two vertical coordinates which define the horizontal stip.
    """
    assert is_k_shape(p, k)
    r = k_row_lengths(p, k)
    # pad with a row of infinite length and a row of length 0
    r = [float('inf')] + r + [0]
    y_min = max([j for j in range(0, len(r)) if r[j] > width])
    y_max = min([j for j in range(0, len(r)) if r[j] < width]) - 1
    return (y_min, y_max)

def v_bounds(p, k, height):
    """ Recall "V_i".  This is the vertical analog of h_bounds. """
    return h_bounds(p.conjugate(), k, height)

def kShape_is_k_reducible_by_rectangle(p, k, (a,b)):
    """ Checks if the k-shape is k-reducible for a k-rectangle of specific dimensions a x b.

    See Proposition 3.8 in Combinatorics of k-shapes and Genocchi
numbers
    """
    assert is_k_shape(p, k)
    assert a + b - 1 == k or a + b - 1 == k - 1
    # get intersection H_a \cap V_b \cap k_rim
    rim = k_rim(p, k)
    (y_min, y_max) = h_bounds(p, k, a)
    (x_min, x_max) = v_bounds(p, k, b)
    intersection_rim = [(x,y) for (x,y) in rim if x_min <= x <= x_max and y_min <= y <= y_max]
    # check condition (iii) of Proposition 3.8
    if not intersection_rim:
        return False
    else:
        # min_y is DIFFERENT than y_min
        min_y = intersection_rim[0][1]
        max_y = intersection_rim[-1][1]
        return max_y - min_y >= b

def kShape_is_k_reducible2(p, k):
    rect_dim_list = k_rectangle_dimension_list(k) + k_rectangle_dimension_list(k-1)
    for (a, b) in rect_dim_list:
        if kShape_is_k_reducible_by_rectangle(p, k, (a,b)):
            return True
    return False

def kShape_is_k_reducible(s, k, method=2):
    if method == 2:
        return kShape_is_k_reducible2(s, k)
    else:
        raise ValueError('Unknown reducibility method.')

def kShape_is_k_irreducible(s, k, method=2):
    """ A k-shape is called __k-irreducible__ if it is not k-reducible. """
    return not kShape_is_k_reducible(s, k, method)
# END k-shape methods.

def get_k_irreducible_k_shapes(k, method=2):
    # The k-row-shape has at most k rows of length 0, k-1 rows of length 1, ..., 0 rows of length k.  And 0 rows of length greater than k.  Hence the k-row-shape has an upper bound of k*(k-1)/2 rows.  The same goes for the k-col-shape.
    bound = (k-1)*k/2
    n_bound = bound**2
    ptns = []
    for n in range(0, n_bound+1):
        ptns += Partitions(n, max_length=bound, max_part=bound)
        k_irr_k_shapes = [p for p in ptns if is_k_shape(p, k) and kShape_is_k_irreducible(p, k, method)]
    return k_irr_k_shapes

