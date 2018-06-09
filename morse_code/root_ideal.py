# -*- coding: utf-8 -*-
"""
Sage does *not* have a builtin 'RootIdeal' object.  *This* module contains a RootIdeal class and useful functions pertaining to root ideals:

REFERENCES:

.. [cat] `Catalan functions and k-schur positivity <https://arxiv.org/abs/1804.03701>`_
.. [scat] Skew-linked Catalan functions and k-schur positivity.  Jonah Blasiak, Jennifer Morse, Anna Pun, and Daniel Summers.  Not to be confused with 'Catalan functions and k-schur positivity.'
"""
from sage.all import *
import partition as P
import skew_partition as SP

# HELPERS
def is_weakly_decreasing(li):
    return all(li[i] >= li[i+1] for i in range(len(li)-1))

def is_strictly_decreasing(li):
    return all(li[i] > li[i+1] for i in range(len(li)-1))

def get_n_from_root_ideal(root_ideal):
    return max(c for (r,c) in root_ideal) + 1

def get_dim(n, ri_list):
    if n is not None:
        return n
    for ri in ri_list:
        # To safeguard against user error, perhaps we should ALWAYS require user to input n.
        if ri and isinstance(ri, RootIdeal):
            return get_n_from_root_ideal(ri)
    raise Exception('There is no way to figure out the size of the staircase that the root ideals fall in.  Please supply n.')


# RootIdeal stuff
class RootIdeal(list):
    """ An upper root ideal.

    Consider the k-1 staircase partition `[k-1, k-2, \\ldots, 1]` positioned in the upper-right corner of a `k` x `k` grid.  The cells in the grid are labeled with (row_index, col_index) 0-based coordinates.  Now consider any right-justified subpartition of the staircase partition.  This is a RootIdeal.  However, it is expressed not as a partition but as a list of the cells it contains.

    For example, the partition `[3, 1]` in the 7 x 7 grid is the root ideal `[(0,4), (0,5), (0,6), (1,6)]`.

    See Definition 2.1 of [cat]_ for more.
    """
    def __hash__(self):
        return hash(tuple(sorted(self)))

def bump_path_piece(sp, start_row_index, blocked_rows=set()):
    # Helper
    # this algo find the correct "L" piece of the path, where the bottom right cell is cell1, the bottom left is cell2, and the top left is cell3
    # Returns (top_row_index, is_end) which are the row index of cell3 and whether or not we 'broke free' out of the top left cell of the skew-partition, respectively.
    col_index2 = SP.left(sp, start_row_index)
    row_index3 = SP.top(sp, col_index2) + 1
    while row_index3 in blocked_rows:
        row_index3 += 1
    # CATTY-CORNER ONLY line:
    max_row_index = len(sp.outer()) - 1
    if row_index3 > max_row_index:
        return None, True
    else:
        return row_index3, False
def bump_path(sp, row_index, blocked_rows=set()):
    # helper
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
    # Given the dimension of the square n and the selected rows, output the root ideal
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
    Given a SkewPartition `sp`, return the removable roots of the corresponding 'min' or 'max' root ideal.

    type: 'min' or 'max'
    returns: A list of removable roots in order.
    """
    def type_shift(x, type):
        if type == 'max':
            return x
        elif type == 'min':
            return x - 1
        else:
            raise ValueError('Bad type.')
    assert SP.is_linked(sp)
    mu = Partition(sp.column_lengths())
    eta = sp.inner()
    rmvble_roots = []
    for i in range(0, len(eta)):
        mu_index = type_shift(eta[i], type)
        rmvble_root = (i, mu[mu_index] + i)
        rmvble_roots.append(rmvble_root)
    return rmvble_roots

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
    """ Given the removable roots `corners` of a root ideal and the size length `n` of the n x n grid, return the root ideal itself. """
    ptn = removable_roots_to_partition(corners, n)
    ri = partition_to_root_ideal(ptn, n)
    return ri

def skew_partition_to_root_ideal(sp, type='max', method='removable roots'):
    """ Given a SkewPartition `sp` and a type of root ideal ('max' or 'min'), return the corresponding root ideal. """
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
    return RootIdeal(root_ideal)

def RootIdeal_next(ri, min=[], max=None, n=None, type='strict'):
    # figure out dimension of square
    n = get_dim(n, [ri, min, max])
    ptn = root_ideal_to_partition(ri)
    min_ptn = root_ideal_to_partition(min)
    max_ptn = root_ideal_to_partition(max)
    if type in ('strict', 'rational'):
        type = 'strictly decreasing'
    next_ptn = P.next(ptn, min=min_ptn, max=max_ptn, type=type)
    next_ri = partition_to_root_ideal(next_ptn, n)
    return next_ri

def skew_partition_to_root_ideals(sp, type='strict'):
    """ Given a skew partition `sp`, find the corresponding set (but given as a list here) of root ideals.
    """
    # We could change this to an iterator if users may not want all the root ideals.
    min_ri = skew_partition_to_root_ideal(sp, type='min')
    max_ri = skew_partition_to_root_ideal(sp, type='max')
    n = len(sp.outer())
    next_func = lambda ri: RootIdeal_next(ri, min=min_ri, max=max_ri, n=n, type=type)
    ris = generate_path(next_func, min_ri)
    return ris

def down(ri, row_index):
    """ Given a root ideal `ri` and a starting position `row_index`, move right on that row until you hit the root ideal (you are now standing ontop of a cell of the root ideal), then move straight down until you hit the diagonal, and return the new index.
    """
    # Note: I am assuming the cells in the root ideal are IN ORDER with y coordinates weakly increasing, and for fixed y, x strictly increasing
    for (r,c) in ri:
        if r == row_index:
            return c
    return None

def up(root_ideal, col_index):
    """ Same as :meth:`down`, but this time you start in the *column* indicated by `column_index`, and move *up* until you hit the root ideal, then move *left* until you hit the diagonal.
    """
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
    """ Given a starting row index `start_index`, perform :meth:`down` operations repeatedly until you can't anymore.  Returns the resulting sequence of indecis as a list.
    """
    next_func = lambda index: down(root_ideal, index)
    return generate_path(next_func, start_index)

def up_path(root_ideal, start_index):
    """ Same as :meth:`down_path`, but uses a *column* index to start with, and applies *up* operations repeatedly.
    """
    next_func = lambda index: up(root_ideal, index)
    return generate_path(next_func, start_index)

def top(root_ideal, start_index):
    """ Given a column index `start_index`, look at it's :meth:`up_path` and return the final index.
    """
    return up_path(root_ideal, start_index)[-1]

def bottom(root_ideal, start_index):
    """ Given a row index `start_index`, look at it's :meth:`down_path` and return the final index.
    """
    return down_path(root_ideal, start_index)[-1]

def down_path_column_lengths_part(root_ideal, ptn, start_index):
    """ This is `\\mu_i` in Definition 2.3 of [scat]_. """
    return sum(ptn[j] for j in down_path(root_ideal, start_index))
def down_path_column_lengths(root_ideal, ptn):
    """ This is the column shape `\\mu'` as defined by Definition 2.3 of [scat]_.  It is also introduced in the second paragraph of the overview as `\\mathfrak{cs}(\\Psi, \\lambda)`. """
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
                mu.append(down_path_column_lengths_part(root_ideal, ptn, index))
                # remove indecis from future draws
                dpath = down_path(root_ideal, index)
                indecis_available -= set(dpath)
    return Partition(mu)

def root_ideal_to_partition(root_ideal):
    """ Given a root ideal (list of cells), return the corresponding partition (the row shape of the root ideal).
    """
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
    """ Given a partition and the size of the square, return the corresponding root ideal.  (This is the inverse function to :meth:`root_ideal_to_partition` in the context of an `n` x `n` grid.)
    """
    if ptn is None or ptn == False:
        return ptn
    root_ideal = []
    for r, part in enumerate(ptn):
        root_ideal += [(r, c) for c in range(n-part, n)]
    return RootIdeal(root_ideal)

def is_strict(ri):
    """ Given a root ideal `ri`, check to see if it is a *strict root ideal*, as defined in Example 2.4 of [scat]_.  This merely means that it's corresponding partition is strictly decreasing!
    """
    ptn = root_ideal_to_partition(ri)
    return is_strictly_decreasing(ptn)

def complement(ri, n=None):
    """ Given a root ideal (could be upper or lower), return it's complement in the upper-staircase-shape, the result being a root ideal.
    """
    n = get_dim(n, [ri])
    p_staircase = Partition(list(range(n-1, 0, -1)))
    ri_staircase = partition_to_root_ideal(p_staircase, n)
    ri_complement_set = set(ri_staircase) - set(ri)
    ri_complement = sorted(ri_complement_set)
    return ri_complement

def partition_to_k_Schur_root_ideal(ptn, k, n=None):
    """ Given a `k`-bounded partition `ptn` and the dimension `n` of the `n` x `n` grid, return the corresponding `k`-Schur root ideal.
    """
    ptn = Partition(ptn)
    if n is None:
        n = len(ptn)
    assert P.is_k_bounded(ptn, k)
    assert len(ptn) <= n
    ri = []
    for i, part in enumerate(ptn):
        ri += [(i,j) for j in range(k - part + i + 1, n)]
    return ri





