# -*- coding: utf-8 -*-
r"""
Sage does *not* have a builtin 'RootIdeal' object.  *This* module contains a RootIdeal class and useful functions pertaining to root ideals:

REFERENCES:

.. [cat] `Catalan functions and k-schur positivity <https://arxiv.org/abs/1804.03701>`_
.. [scat] Skew-linked Catalan functions and k-schur positivity.  Jonah Blasiak, Jennifer Morse, Anna Pun, and Daniel Summers.  Not to be confused with 'Catalan functions and k-schur positivity.'
"""
from sage.all import *
import partition
import skew_partition
# ^*^ sphinx insert ^*^

# HELPERS
def is_weakly_decreasing(li):
    return all(li[i] >= li[i+1] for i in range(len(li)-1))

def is_strictly_decreasing(li):
    return all(li[i] > li[i+1] for i in range(len(li)-1))

def generate_path(next_func, start):
    path = [start]
    while True:
        next_ = next_func(path[-1])
        if next_ is not None:
            path.append(next_)
        else:
            break
    return path

def staircase_shape(n):
    r""" Given `n`, return the composition `[n-1, n-2, \ldots, 0]` commonly denoted `\rho`.

    Yes, this INCLUDES a 0 at the end!
    """
    return Composition(range(n - 1, -1, -1))

def bump_path_piece(sp, start_row_index, blocked_rows=set()):
    # Helper
    # this algo find the correct "L" piece of the path, where the bottom right cell is cell1, the bottom left is cell2, and the top left is cell3
    # Returns (top_row_index, is_end) which are the row index of cell3 and whether or not we 'broke free' out of the top left cell of the skew-partition, respectively.
    col_index2 = skew_partition.left(sp, start_row_index)
    row_index3 = skew_partition.top(sp, col_index2) + 1
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

def selected_rows_to_maximum_root_ideal(n, selected_indices):
    # Given the dimension of the square n and the selected rows, output the root ideal
    root_ideal_cells = []
    selected_indices = set(selected_indices)
    permitted_col_indices = set(range(n)) - selected_indices
    for i in range(n):
        if i in selected_indices:
            if permitted_col_indices:
                smallest_unblocked_index = min(permitted_col_indices)
                root_ideal_cells += [(i, j) for j in range(smallest_unblocked_index, n)]
                permitted_col_indices.remove(smallest_unblocked_index)
                selected_indices.add(smallest_unblocked_index)
    return RootIdeal(root_ideal_cells)

def skew_partition_to_removable_roots(sp, type='max'):
    r"""
    Given a SkewPartition ``sp``, return the removable roots of the corresponding 'min' or 'max' root ideal.

    INPUTS:

    - ``sp`` -- a SkewPartition

    OPTIONAL INPUTS:

    - ``type`` -- (default ``'max'``) the type of root ideal you want to use.  ``'min'`` is the minimum root ideal (as far as containment goes) and ``'max'`` is the maximum root ideal.

    OUTPUT:

    A list of removable roots in order.

    """
    def type_shift(x, type):
        if type == 'max':
            return x
        elif type == 'min':
            return x - 1
        else:
            raise ValueError('Bad type.')
    assert skew_partition.is_linked(sp)
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

# RootIdeal stuff
def is_roots(obj):
    # Dirty indicator of whether object is roots (is it an iterable of pairs of natural numbers).
    try:
        iter(obj)
    except:
        return False
    if not all(isinstance(el, tuple) and len(el) == 2 and el[0] in NonNegativeIntegerSemiring() and el[1] in NonNegativeIntegerSemiring() for el in obj):
        return False
    return True


class RootIdeal(list):
    r""" An upper root ideal.

    Consider the k-1 staircase partition `[k-1, k-2, \ldots, 1]` positioned in the upper-right corner of a `k` x `k` grid.  The cells in the grid are labeled with (row_index, col_index) 0-based coordinates.  Now consider any right-justified subpartition of the staircase partition.  This is a RootIdeal.  However, it is expressed not as a partition but as a list of the cells it contains.

    For example, the partition `[3, 1]` in the 7 x 7 grid is the root ideal `[(0,4), (0,5), (0,6), (1,6)]`.

    See Definition 2.1 of [cat]_ for more.

    EXAMPLES::

        sage: ri = RootIdeal([(0,4), (0,5), (0,6), (1,6)])

    """
    def __init__(self, lis, n=None):
        # validate the roots
        assert is_roots(lis)
        # figure out n
        if n is not None:
            self.n = n
        elif lis:
            self.n = max(c for (r,c) in lis) + 1
        else:
            self.n = 1
        # normalize the roots
        lis = sorted(lis)
        list.__init__(self, lis)

    def __hash__(self):
        return hash(tuple(self))

    def next(self, min=[], max=None, type='strict'):
        # figure out dimension of square
        n = self.n
        ptn = self.to_partition()
        min_ptn = RootIdeal(min, n=n).to_partition()
        max_ptn = RootIdeal(max, n=n).to_partition() if max is not None else None
        if type in ('strict', 'rational'):
            type = 'strictly decreasing'
        next_ptn = partition.next(ptn, min=min_ptn, max=max_ptn, type=type)
        next_ri = RootIdeals().init_from_partition(next_ptn, n)
        return next_ri

    def down(ri, row_index):
        r""" Given a root ideal `ri` and a starting position 'row_index', move right on that row until you hit the root ideal (you are now standing ontop of a cell of the root ideal), then move straight down until you hit the diagonal, and return the new index.

        The picture below represents the root ideal used in the example.

        .. image:: _static/root-ideal.JPG
            :width: 180px
            :align: center
            :alt: The root ideal [(0,2), (0,3), (0,4), (1,3), (1,4), (2,4)]

        EXAMPLES::

            sage: ri = RootIdeals().init_from_partition([3, 2, 1], 5)
            sage: ri.down(0)
            2
            sage: ri.down(1)
            3
            sage: ri.down(2)
            4
            sage: ri.down(3) == None
            True
            sage: ri.down(4) == None
            True

        """
        # Note: I am assuming the cells in the root ideal are IN ORDER with y coordinates weakly increasing, and for fixed y, x strictly increasing
        for (r,c) in ri:
            if r == row_index:
                return c
        return None

    def up(root_ideal, col_index):
        r""" Same as :meth:`down`, but this time you start in the *column* indicated by 'column_index', and move *up* until you hit the root ideal, then move *left* until you hit the diagonal.

        The picture below represents the root ideal used in the example.

        .. image:: _static/root-ideal.JPG
            :width: 180px
            :align: center
            :alt: The root ideal [(0,2), (0,3), (0,4), (1,3), (1,4), (2,4)]

        EXAMPLES::

            sage: ri = RootIdeals().init_from_partition([3, 2, 1], 5)
            sage: ri.up(0) == None
            True
            sage: ri.up(1) == None
            True
            sage: ri.up(2)
            0
            sage: ri.up(3)
            1
            sage: ri.up(4)
            2

        """
        for (r,c) in reversed(root_ideal):
            if c == index:
                return r
        return None

    def down_path(root_ideal, start_index):
        r""" Given a starting row index 'start_index', perform :meth:`down` operations repeatedly until you can't anymore.  Returns the resulting sequence of indices as a list.  (See [cat]_ Definition 5.2 for more)

        The picture below represents the root ideal used in the example, and the path drawn on the picture depicts the down path for ``start_index`` 0 specifically.

        .. image:: _static/bottom.JPG
            :width: 180px
            :align: center
            :alt: The root ideal [(0,2), (0,3), (0,4), (1,3), (1,4), (2,4)]

        EXAMPLES::

            sage: ri = RootIdeals().init_from_partition([3, 2, 1], 5)
            sage: ri.down_path(0)
            [0, 2, 4]
            sage: ri.down_path(1)
            [1, 3]
            sage: ri.down_path(2)
            [2, 4]
            sage: ri.down_path(3)
            [3]
            sage: ri.down_path(4)
            [4]

        """
        next_func = lambda index: root_ideal.down(index)
        return generate_path(next_func, start_index)

    def up_path(root_ideal, start_index):
        r""" Same as :meth:`down_path`, but uses a *column* index to start with, and applies *up* operations repeatedly.

        The picture below represents the root ideal used in the example.

        .. image:: _static/root-ideal.JPG
            :width: 180px
            :align: center
            :alt: The root ideal [(0,2), (0,3), (0,4), (1,3), (1,4), (2,4)]

        EXAMPLES::

            sage: ri = RootIdeals().init_from_partition([3, 2, 1], 5)
            sage: ri.up_path(0)
            [0]
            sage: ri.up_path(1)
            [1]
            sage: ri.up_path(2)
            [2, 0]
            sage: ri.up_path(3)
            [3, 1]
            sage: ri.up_path(4)
            [4, 2, 0]

        """
        next_func = lambda index: root_ideal.up(index)
        return generate_path(next_func, start_index)

    def top(root_ideal, start_index):
        r""" Given a column index 'start_index', look at it's :meth:`up_path` and return the final index.

        The picture below represents the root ideal used in the example.

        .. image:: _static/root-ideal.JPG
            :width: 180px
            :align: center
            :alt: The root ideal [(0,2), (0,3), (0,4), (1,3), (1,4), (2,4)]

        EXAMPLES::

            sage: ri = RootIdeals().init_from_partition([3, 2, 1], 5)
            sage: ri.top(0)
            0
            sage: ri.top(1)
            1
            sage: ri.top(2)
            0
            sage: ri.top(3)
            1
            sage: ri.top(4)
            0

        """
        return root_ideal.up_path(start_index)[-1]

    def bottom(root_ideal, start_index):
        r""" Given a row index 'start_index', look at it's :meth:`down_path` and return the final index.

        The picture below represents the root ideal used in the examples, and the path drawn on the picture depicts the down path for index 0 specifically, demonstrating that ``bottom(0)`` should be 4.

        .. image:: _static/bottom.JPG
            :width: 180px
            :align: center
            :alt: The root ideal [(0,2), (0,3), (0,4), (1,3), (1,4), (2,4)]

        EXAMPLES::

            sage: ri = RootIdeals().init_from_partition([3, 2, 1], 5)
            sage: ri.bottom(0)
            4
            sage: ri.bottom(1)
            3
            sage: ri.bottom(2)
            4
            sage: ri.bottom(3)
            3
            sage: ri.bottom(4)
            4

        """
        return root_ideal.down_path(start_index)[-1]

    def down_path_column_lengths_part(root_ideal, ptn, start_index):
        r""" This is `\mu_i` in Definition 2.3 of [scat]_.

        This exists mainly as a helper function for :meth:`down_path_column_lengths`.
        """
        return sum(ptn[j] for j in root_ideal.down_path(start_index))
    def down_path_column_lengths(self, ptn):
        r""" This is the column shape `\mu'` as defined by Definition 2.3 of [scat]_.  It is also introduced in the second paragraph of the overview as `\mathfrak{cs}(\Psi, \lambda)`.

        In Example 2.4 of [scat]_, the following

        ..  image:: _static/example2.4.png
            :align: center
            :alt: The root ideal [(0,1), (0,2), (0,3), (0,4), (0,5), (1,4), (1,5), (2,4), (2,5), (3,4), (3,5)] and the partition 7 6 5 2 2 2

        depicts the root ideal in red and the partition on the diagonal.

        EXAMPLES::

            sage: ri = RootIdeals().init_from_partition([5, 2, 2, 2], 6)
            sage: ptn = [7, 6, 5, 2, 2, 2]
            sage: ri.down_path_column_lengths(ptn)
            [15, 7, 4, 2]

        This is also the lengths of the bounce paths in [cat]_ Definition 5.2.

    """
        if not self:
            mu = ptn
        else:
            mu = []
            # n is the side length of the square
            n = self.n
            indices_available = set(range(0, n))
            for index in range(0, n):
                if index in indices_available:
                    # add the kthing to mu
                    mu.append(self.down_path_column_lengths_part(ptn, index))
                    # remove indices from future draws
                    dpath = self.down_path(index)
                    indices_available -= set(dpath)
        return Partition(mu)

    def to_partition(root_ideal):
        r""" Given a root ideal (list of cells), return the corresponding partition (the row shape of the root ideal).

        The red part of the following picture (please ignore the diagonal) represents the root ideal [(0, 1), (0, 2), (0, 3), (0, 4), (0, 5), (1, 4), (1, 5), (2, 4), (2, 5), (3, 4), (3, 5)].

        .. image:: _static/Ksi.png
            :align: center
            :alt: The root ideal [(0, 1), (0, 2), (0, 3), (0, 4), (0, 5), (1, 4), (1, 5), (2, 4), (2, 5), (3, 4), (3, 5)] inside a 6 by 6 grid.

        But it can also be interpreted as the partition 5 2 2 2 (in the Hebrew convention).  Therefore, ``to_partition()`` acting on the root ideal will output 5 2 2 2.

        EXAMPLES::

            sage: ri = RootIdeal([(0, 1), (0, 2), (0, 3), (0, 4), (0, 5), (1, 4), (1, 5), (2, 4), (2, 5), (3, 4), (3, 5)])
            sage: ri.to_partition()
            [5, 2, 2, 2]

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

    def is_strict(ri):
        r""" Given a root ideal ``ri``, check to see if it is a *strict root ideal*, as defined in Example 2.4 of [scat]_.  This merely means that it's corresponding partition is strictly decreasing!

        In the following images, ignore the index on the diagonal and look only at the root ideal in red.

        .. image:: _static/Phi.png
            :align: center
            :alt: The root ideal corresponding to the partition  8 6 5 3 2

        EXAMPLES::

            sage: ri = RootIdeals().init_from_partition([8, 6, 5, 3, 2], 9)
            sage: ri.is_strict()
            True

        .. image:: _static/Ksi.png
            :align: center
            :alt: The root ideal corresponding to the partition  5 2 2 2

        EXAMPLES::

            sage: ri = RootIdeals().init_from_partition([5, 2, 2, 2], 6)
            sage: ri.is_strict()
            False

        """
        ptn = ri.to_partition()
        return is_strictly_decreasing(ptn)

    def complement(self):
        r""" Given a root ideal (or just an iterable of roots), return it's complement in the upper-staircase-shape, the result being a root ideal (or just an iterable of roots).

        INPUTS:

        - ``ri`` -- a root ideal

        OPTIONAL INPUTS:

        - ``n`` -- (default ``None``) the side length of the n x n box you want the complement to be taken over.

        For example, the two root ideals depicted below are complements of each other:

        .. image:: _static/root-ideal.JPG
            :width: 180px
            :align: center
            :alt: The root ideal [(0,2), (0,3), (0,4), (1,3), (1,4), (2,4)]

        .. image:: _static/root-ideal-complement.JPG
            :width: 180px
            :align: center
            :alt: The root ideal [(0,1), (1,2), (2,3), (3,4)]

        EXAMPLES::

            sage: ri1 = RootIdeal([(0,2), (0,3), (0,4), (1,3), (1,4), (2,4)])
            sage: ri2 = RootIdeal([(0,1), (1,2), (2,3), (3,4)])
            sage: ri1.complement() == ri2
            True
            sage: ri2.complement() == ri1
            True

        """
        n = self.n
        ri_staircase = RootIdeals().init_staircase(n)
        ri_complement_set = set(ri_staircase) - set(self)
        ri_complement = RootIdeal(ri_complement_set)
        return ri_complement


class RootIdeals:
    def init_from_removable_roots(self, corners, n):
        r""" Given the removable roots ``corners`` of a root ideal and the size length `n` of the `n` x `n` grid, return the root ideal itself.

        For example, the root ideal [(0, 1), (0, 2), (0, 3), (0, 4), (0, 5), (1, 4), (1, 5), (2, 4), (2, 5), (3, 4), (3, 5)] in the `6` x `6` grid shown below in red (please ignore the diagonal) has removable roots `(0, 1)` and `(3, 4)`.

        .. image:: _static/Ksi.png
            :align: center
            :alt: The root ideal corresponding to the partition  5 2 2 2

        EXAMPLES::

            sage: removable_roots_to_root_ideal({(0, 1), (3, 4)}, 6)
            [(0, 1), (0, 2), (0, 3), (0, 4), (0, 5), (1, 4), (1, 5), (2, 4), (2, 5), (3, 4), (3, 5)]

        """
        ptn = removable_roots_to_partition(corners, n)
        ri = self.init_from_partition(ptn, n)
        return ri

    def init_from_skew_partition(self, sp, type='max', method='removable roots'):
        r""" Given a SkewPartition `sp` and a type of root ideal ('max' or 'min'), return the corresponding root ideal.

        A type of ``'min'`` returns `\Phi(\lambda, \mu)` while a type of ``'max'`` returns `\Phi^+(\lambda, \mu)` as notated in [scat]_ at the bottom of page 1.

        """
        if method == 'removable roots':
            corners = skew_partition_to_removable_roots(sp, type)
            n = len(sp.outer())
            root_ideal = self.init_from_removable_roots(corners, n)
        elif method == 'bounce':
            if type != 'max':
                raise Exception('The bounce method can only yield the maximum root ideal (type=max).')
            selected_indices = skew_partition_to_selected_rows(sp)
            n = len(sp.outer())
            root_ideal = selected_rows_to_maximum_root_ideal(n, selected_indices)
        else:
            raise ValueError('Unknown method.')
        return RootIdeal(root_ideal)

    def init_all_from_skew_partition(self, sp, type='strict'):
        r""" Given a skew partition `sp`, find the corresponding set (but given as a list here) of root ideals.

        (This is the set `\{\Psi \in \Delta^+(\mathfrak{R}) \mid \Phi(\lambda, \mu) \subset \Psi \subset \Phi^+(\lambda, \mu)\} = [(\lambda, \mu)]` found in [scat]_ at the bottom of page 1.)

        """
        # We could change this to an iterator if users may not want all the root ideals.
        min_ri = self.init_from_skew_partition(sp, type='min')
        max_ri = self.init_from_skew_partition(sp, type='max')
        n = len(sp.outer())
        next_func = lambda ri: ri.next(min=min_ri, max=max_ri, type=type)
        ris = generate_path(next_func, min_ri)
        return ris

    def init_k_schur_from_partition(self, ptn, k, n=None):
        r""" Given a `k`-bounded partition `ptn = \mu` and the dimension `n` of the `n` x `n` grid, return the corresponding `k`-Schur root ideal `\Delta^k(\mu)`, as defined in [cat]_ Definition 2.2 as

        .. math::

            \Delta^k(\mu) = \{(i,j) \in \Delta^+_n \mid k - \mu_i + i < j \}

        The following diagram depicts the `k`-bounded partition on the diagonal and the resulting `k`-schur root ideal.

        .. image:: _static/k-schur-root-ideal.JPG
            :height: 200px
            :align: center
            :alt: The partition 3 3 2 2 2 and the root ideal [(0, 2), (0, 3), (0, 4), (0, 5), (0, 6), (0, 7), (1, 3), (1, 4), (1, 5), (1, 6), (1, 7), (2, 5), (2, 6), (2, 7), (3, 6), (3, 7), (4, 7)]

        EXAMPLES::

            sage: k_ri = RootIdeals().init_k_schur_from_partition([3, 3, 2, 2, 2], 4, n=8)
            sage: k_ri.to_partition()
            [6, 5, 3, 2, 1]

        """
        ptn = Partition(ptn)
        if n is None:
            n = len(ptn)
        assert partition.is_k_bounded(ptn, k)
        assert len(ptn) <= n
        ri = []
        for i, part in enumerate(ptn):
            ri += [(i,j) for j in range(k - part + i + 1, n)]
        return RootIdeal(ri)

    def init_from_partition(self, ptn, n):
        r""" Given a partition and the size of the square, return the corresponding root ideal.  (This is the inverse function to :meth:`RootIdeal.to_partition` in the context of an `n` x `n` grid.)

        The red part of the following picture (please ignore the diagonal) can be interpreted as the partition 5 2 2 2 (in the Hebrew convention):

        .. image:: _static/Ksi.png
            :align: center
            :alt: The root ideal [(0, 1), (0, 2), (0, 3), (0, 4), (0, 5), (1, 4), (1, 5), (2, 4), (2, 5), (3, 4), (3, 5)] inside a 6 by 6 grid.

        Therefore the partition 5 2 2 2 with `n=6` corresponds to the root ideal [(0, 1), (0, 2), (0, 3), (0, 4), (0, 5), (1, 4), (1, 5), (2, 4), (2, 5), (3, 4), (3, 5)].

        EXAMPLES::

            sage: RootIdeals().init_from_partition([5, 2, 2, 2], 6)
            [(0, 1), (0, 2), (0, 3), (0, 4), (0, 5), (1, 4), (1, 5), (2, 4), (2, 5), (3, 4), (3, 5)]

        """
        if ptn is None or ptn == False:
            return ptn
        ptn = Partition(ptn)
        root_ideal = []
        for r, part in enumerate(ptn):
            root_ideal += [(r, c) for c in range(n-part, n)]
        return RootIdeal(root_ideal)

    def init_staircase(self, n):
        r""" Given `n`, return the root ideal commonly denoted `\Delta^+`, which is the maximum possible root ideal in an `n` x `n` grid.

        EXAMPLES::

            sage: RootIdeals().init_staircase(3)
            [(0,1), (0,2), (1,2)]
            sage: RootIdeals().init_staircase(4)
            [(0,1), (0,2), (0,3), (1,2), (1,3), (2,3)]

        """
        return self.init_from_partition(staircase_shape(n), n)
