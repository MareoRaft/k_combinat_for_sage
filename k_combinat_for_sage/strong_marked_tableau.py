# -*- coding: utf-8 -*-
r"""
Sage has have a builtin `StrongTableau <https://doc.sagemath.org/html/en/reference/combinat/sage/combinat/k_tableau.html>`_ object.  *This* module contains useful functions pertaining to Standard Marked Tableau (SMT), the most standard form of a StrongTableau (a strong tableau with standard weight).

AUTHORS:

- Matthew Lancellotti (2018): Initial version

REFERENCES:
"""

#*****************************************************************************
#  Copyright (C) 2018 Matthew Lancellotti <mvlancellotti@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.all import *

from partition import *
import skew_partition
# ^*^ sphinx insert ^*^


def shape_cell_indices(shape):
    # r"""
    # Takes a shape (partition or list of lists) and returns the
    # indices of the cells in the shape.
    # """
    shape = Partition(shape)
    return shape.cells()


def tableau_cell_indices(tab):
    # r"""
    # Takes a tableau (or list of lists) and returns the indices
    # of the cells in the tableau.
    # """
    tab = Tableau(tab)
    return shape_cell_indices(tab.shape())


def tableau_contains(outer, inner):
    # r"""
    # Checks to see if `inner` tableau is contained in the `outer`
    # tableau exactly. Somewhat dangerous function.
    # """
    ol = outer.to_list()
    il = inner.to_list()
    result = True
    indices = tableau_cell_indices(inner)
    for (row_index, col_index) in indices:
        if il[row_index][col_index] != ol[row_index][col_index] and il[row_index][col_index] is not None:
            result = False
            break
    return result


def strong_tableau_quotient2(outer_tab, inner_shape):
    # r"""
    # Takes a strong tableau and skews it; that is, it sets all the
    # cells in `outer_tab` that are contained in `inner_shape` to `None`
    # and returns a `Tableau` object.
    # """
    ol = outer_tab.to_list()
    indices = shape_cell_indices(inner_shape)
    for (row_index, col_index) in indices:
        ol[row_index][col_index] = None
    return Tableau(ol)


def strong_tableau_quotient(outer, inner):
    # r"""
    # Takes a strong tableau and skews it; that is, it sets all the
    # cells in `outer_tab` that are contained in `inner_shape` to `None`
    # and returns a `StrongTableau` object.
    # """
    assert tableau_contains(outer, inner)
    ol = outer.to_list()
    indices = tableau_cell_indices(inner)
    for (row_index, col_index) in indices:
        ol[row_index][col_index] = None
    return StrongTableau(ol, outer.k)


def strong_tableau_has_row_marking(tab, row_index):
    r""" Returns whether ``tab`` has a row marking in row ``row_index``.

    Checks if a strong tableau ``tab`` has a row marking in row
    ``row_index`` (row indexing starts at 0) and returns ``True``
    if row_marking is present; ``False`` otherwise.

    ..  SEEALSO::

        :meth:`is_row_markable`
    """
    # WARNING: Indices start at 0
    if len(tab) <= row_index:
        return False
    row = tab[row_index]
    result = (True in [i is not None and i < 0 for i in row])
    return result


def add_skew_tabs(tab1, tab2):
    # Delete this
    # Dangerous if overlapping entries / has 0's
    m1 = matrix(RationalField(), tab1.to_list())
    m2 = matrix(RationalField(), tab2.to_list())
    added = m1+m2
    tab_list = [map(lambda i: None if i == 0 else i, row) for row in added]
    return Tableau(tab_list)


def add_skew_tab_to_tab(tab, skew_tab):
    # r"""
    # Given a tableau `tab` as a base, this function layers a skew tableau
    # `skew_tab` on top of it.
    # """
    skew_list = skew_tab.to_list()
    # Necessary because skew strong tableaux does not get item correctly?
    tab_list = tab.to_list()
    for (i, j) in tableau_cell_indices(tab):
        skew_list[i][j] = tab_list[i][j]
    return Tableau(skew_list)


def shape_to_homogeneous_tab(outer_shape, entry, inner_shape=[]):
    li = []
    for i in range(len(outer_shape)):
        row = []
        if len(inner_shape) > i:
            row = row + [None]*inner_shape[i] + \
                [entry]*(outer_shape[i]-inner_shape[i])
        else:
            row = row + [entry]*outer_shape[i]
        li.append(row)
    return Tableau(li)


def std_strong_tab_from_core_sequence(core_sequence, k, markings):
    # 'markings' are (row,col) coordinates
    if len(core_sequence) == 0:
        return StrongTableaux(k, [], []).an_element()
    wt = [1]*(len(core_sequence)-1)
    core_zero = core_sequence[0]
    base = Tableau([[None]*i for i in core_zero])
    for core_index in range(1, len(core_sequence)):
        skew = shape_to_homogeneous_tab(
            core_sequence[core_index], core_index, core_sequence[core_index-1])
        base = add_skew_tab_to_tab(base, skew)
    return StrongTableaux.add_marking(base, markings, k, wt)


def _strong_marked_tableau(lis, k):
    # helper to create a strong marked tableau
    st = StrongTableau(lis, k)
    for w in st.weight():
        if w != 1:
            raise ValueError('Not a standard weight SMT.')
    if st.to_list() != st.to_standard_list():
        raise ValueError('Not a standard SMT.')
    return st


def k_coverees1(root, k):
    # THIS FUNCTIONALITY IS ALREADY BUILTIN.  See algorithm 2 of k_coverees.
    # one way to get the k coverees
    root = Core(root, k+1)
    root = root.to_partition()
    # set of coveree candidates (a superset of the coverees)
    candidates = set()

    def add_to_candidates(ptn, is_root=False):
        # add to dictionary if not already there
        if ptn not in candidates:
            candidates.add(ptn)
            # if it is not a leaf, recurse
            if is_root or not ptn.is_core(k+1):
                for (i, j) in ptn.removable_cells():
                    sub_ptn = ptn.remove_cell(i, j)
                    add_to_candidates(sub_ptn)
    add_to_candidates(root, is_root=True)
    # Now that all k+1-cores (and some other things) have been populated in candidates, filter to the ones that really are k+1-cores and have correct k-boundary size.
    coverees = set(ptn for ptn in candidates if ptn.is_core(k+1)
                   and k_size(ptn, k) == k_size(root, k) - 1)
    return coverees


def k_coverees(core, k, algorithm=1):
    # THIS FUNCTIONALITY IS ALREADY BUILTIN.  See algorithm 2 below.
    r""" Given a `k+1`-core, find all sub-`k+1`-cores that have `k`-boundary 1 less than the given. """
    if algorithm == 1:
        return k_coverees1(core, k)
    elif algorithm == 2:
        core = Core(core, k+1)
        coveree_core_list = core.strong_down_list()
        coverees = set(c.to_partition() for c in coveree_core_list)
        return coverees
    else:
        raise ValueError('Unknown algorithm.')


def __go_to_ribbon_head(cells, start_cell):
    # Given the cells of a ribbon or multiple disconnected ribbons, and a starting point, find the head of the ribbon
    if start_cell not in cells:
        raise ValueError('Starting position is not in list of cells.')
    cell = start_cell
    while True:
        right_cell = (cell[0], cell[1] + 1)
        if right_cell in cells:
            cell = right_cell
            continue
        down_cell = (cell[0] - 1, cell[1])
        if down_cell in cells:
            cell = down_cell
            continue
        break
    # nothing to right or down, so we must have reached head of ribbon
    head = cell
    return head


def row_marking_to_marking(outer_core, inner_core, row_marking):
    r""" Convert a row marking to a marking.

    Given the skew shape defined by ``outer_core`` and ``inner_core``, and the row index ``row_marking`` where a marking exists, return the marking `(\text{marking row index}, \text{marking column index})`.

    Note that `\text{marking row index}` mentioned above is simply ``row_marking``.  Therefore, the real usefulness of this function is that it finds the column index of the marking.
    """
    sp = SkewPartition([outer_core, inner_core])
    cells = sp.cells()
    row_indices = set(cell[0] for cell in cells)
    if row_marking not in row_indices:
        raise ValueError('no such row marking')
    start_cell = (row_marking, skew_partition.right(sp, row_marking))
    head = __go_to_ribbon_head(cells, start_cell)
    if head[0] == row_marking:
        return head
    else:
        raise ValueError('no such row marking')


def row_markings_to_markings(core_sequence, row_markings):
    r""" Given a ``core_sequence`` and corresponding ``row_markings`` for each cover of the sequence, convert the row markings to markings and return them.

    Each row marking in ``row_markings`` is merely the row index of where the marking occurs.  The purpose of this function is to convert each row marking to a "marking" which includes the column index.
    """
    assert len(core_sequence) == len(row_markings) + 1
    markings = []
    for index in range(len(row_markings)):
        # get vars
        row_marking = row_markings[index]
        inner_core = core_sequence[index]
        outer_core = core_sequence[index + 1]
        # convert to marking
        marking = row_marking_to_marking(outer_core, inner_core, row_marking)
        markings.append(marking)
    return markings


def is_row_markable(outer_core, inner_core, row_marking):
    r""" Given two cores (typically consecutive cores in a core sequence), see if ``row_marking`` is a possible row_marking of outer_core/inner_core.

    ..  SEEALSO::

        :meth:`strong_tableau_has_row_marking`
    """
    try:
        row_marking_to_marking(outer_core, inner_core, row_marking)
        return True
    except ValueError:
        return False


def k_marked_coverees(core, k, row_marking):
    r""" Given a k+1-core, find all sub-k+1-cores that have k-boundary 1 less than the given with the given row_marking. """
    coverees = k_coverees(core, k)
    marked_coverees = [c for c in coverees
                       if is_row_markable(core, c, row_marking)]
    return set(marked_coverees)


def end_core_to_marked_core_sequences(end_core, k, row_markings):
    r"""
    Return the set of core sequences marked by ``row_markings`` ending in ``end_core``.

    INPUTS:

    - ``end_core`` -- a `k+1`-core

    - ``k`` -- All of the cores in the core sequences are `k+1`-cores.

    - ``row_markings`` -- vector of row-indices indicating the rows of the markings for each core sequence.  (Note that "markings" are row-col-coordinates, while "row_markings" are merely row-coordinates.)

    OUTPUTS:

    - A set of all possible core sequences that end in ``end_core`` and can be marked by the row marking vector ``row_markings``.

    EXAMPLES::

        sage: end_core_to_marked_core_sequences([5, 3, 1], 2, [0, 1, 2, 0, 1])
        {([], [1], [1, 1], [2, 1, 1], [3, 1, 1], [5, 3, 1])}

        sage: end_core_to_marked_core_sequences([5, 3, 1], 2, [1])
        {([3, 1, 1], [5, 3, 1]), ([4, 2], [5, 3, 1])}

    ..  SEEALSO::

        :meth:`end_core_to_strong_marked_tableaux`
    """
    # check inputs
    k = NonNegativeIntegerSemiring()(k)
    end_core = Core(end_core, k+1)
    end_core = end_core.to_partition()
    for row_marking in row_markings:
        NonNegativeIntegerSemiring()(row_marking)
    # find marked core sequences
    sequences = []
    if not row_markings:
        # base case
        sequences.append(tuple([end_core]))
    else:
        # inductive step
        coverees = k_marked_coverees(end_core, k, row_markings[-1])
        for coveree in coverees:
            prefix_sequences = end_core_to_marked_core_sequences(
                coveree, k, row_markings[:-1])
            for prefix_sequence in prefix_sequences:
                sequence = prefix_sequence + tuple([end_core])
                sequences.append(sequence)
    return set(sequences)


def end_core_to_strong_marked_tableaux(end_core, k, row_markings):
    r"""
    Return the set of strong marked tableaux marked by ``row_markings`` ending in ``end_core``.

    INPUTS:

    - ``end_core`` -- a `k+1`-core

    - ``k`` -- All of the cores in the core sequences are `k+1`-cores.

    - ``row_markings`` -- vector of row-indices indicating the rows of the markings for each core sequence.  (Note that "markings" are row-col-coordinates, while "row_markings" are merely row-coordinates.)

    OUTPUTS:

    - A set of all possible strong marked tableau that end in ``end_core`` and can be marked by the row marking vector ``row_markings``.

    EXAMPLES::

        sage: end_core_to_strong_marked_tableaux([5, 3, 1], 2, [0, 1, 2, 0, 1])
        {[[-1, 3, -4, 5, 5], [-2, 5, -5], [-3]]}

        sage: end_core_to_strong_marked_tableaux([5, 3, 1], 2, [1])
        {[[None, None, None, 1, 1], [None, 1, -1], [None]],
         [[None, None, None, None, 1], [None, None, -1], [1]]}

    ..  SEEALSO::

        :meth:`end_core_to_marked_core_sequences`
    """
    core_sequences = end_core_to_marked_core_sequences(
        end_core, k, row_markings)
    smts = set()
    for core_sequence in core_sequences:
        markings = row_markings_to_markings(core_sequence, row_markings)
        smt = std_strong_tab_from_core_sequence(core_sequence, k, markings)
        smts.add(smt)
    return smts
