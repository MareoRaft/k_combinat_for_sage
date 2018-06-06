# -*- coding: utf-8 -*-
"""
Sage does *not* have a builtin 'kShape' object.  *This* module contains useful functions pertaining to `k`-shapes:

REFERENCES:

.. [genocchi] `Combinatorics of k-shapes and Genocchi numbers <https://www.lri.fr/~hivert/PAPER/kshapes.pdf>`_, in FPSAC 2011, Reykjav´k, Iceland DMTCS proc. AO, 2011, 493-504.

"""
from sage.all import *
from partition import *
import skew_partition as SP

# HELPERS
def k_rectangle_dimension_list(k):
    return [(k-i+1, i) for i in range(1, k+1)]

# kShape verifier
def is_k_shape(ptn, k):
    """ A partition is a `k`-*shape* if its `k`-boundary has row-shape and col-shape that are partitions themselves. (Definition 2.1 of [genocchi]_)

    Given a partition `ptn` and a natural number `k`, returns True iff `ptn` is a `k`-shape.

    Given a partition `ptn` *only*, returns True iff there exists some `k \\in [1, n-1]` such that `ptn` is a `k`-shape.

    EXAMPLES::

        sage: is_k_shape(Partition([3, 1]), 1)
        False
        sage: is_k_shape(Partition([3, 1]), 2)
        True

    """
    if k is None:
        # see if it's a k-shape for any k in [1, n-1].
        # (note that every partition is a 0-shape and an n-shape)
        n = ptn.size()
        lis = [is_k_shape(ptn, kk) for kk in range(1, n)]
        return any(lis)
    else:
        k_bdy = ptn.k_boundary(k)
        return SP.is_linked(k_bdy)

#kShape stuff
def h_bounds(p, k, width):
    """ Recall the `H_i` as defined in Definition 3.3 of [genocchi]_.

    Given a natural number `k` (used for the `k`-shape or `k`-boundary) and a width `width`, returns `(y_\\text{min}, y_\\text{max})`, the two vertical coordinates which define the horizontal strip.

    EXAMPLE::

        sage: h_bounds(Partition([10, 7, 4, 2, 2, 2, 1, 1, 1, 1]), k=4, width=3)
        (0, 2)
        sage: h_bounds(Partition([10, 7, 4, 2, 2, 2, 1, 1, 1, 1]), k=4, width=2)
        (2, 3)
        sage: h_bounds(Partition([10, 7, 4, 2, 2, 2, 1, 1, 1, 1]), k=4, width=1)
        (3, 10)

    """
    assert is_k_shape(p, k)
    r = k_row_lengths(p, k)
    # pad with a row of infinite length and a row of length 0
    r = [float('inf')] + r + [0]
    y_min = max([j for j in range(0, len(r)) if r[j] > width])
    y_max = min([j for j in range(0, len(r)) if r[j] < width]) - 1
    return (y_min, y_max)

def v_bounds(p, k, height):
    """ This is `V_i`, the vertical analog of :meth:`h_bounds`.

    EXAMPLE::

        sage: h_bounds(Partition([10, 7, 4, 2, 2, 2, 1, 1, 1, 1]), k=4, width=4)
        (0, 1)
        sage: h_bounds(Partition([10, 7, 4, 2, 2, 2, 1, 1, 1, 1]), k=4, width=3)
        (1, 2)
        sage: h_bounds(Partition([10, 7, 4, 2, 2, 2, 1, 1, 1, 1]), k=4, width=2)
        (2, 2)
        sage: h_bounds(Partition([10, 7, 4, 2, 2, 2, 1, 1, 1, 1]), k=4, width=1)
        (2, 10)

    """
    return h_bounds(p.conjugate(), k, height)

def is_k_reducible_by_rectangle(p, k, (a,b)):
    # Checks if the k-shape is k-reducible for a k-rectangle of specific dimensions a x b.
    # See Proposition 3.8 in Combinatorics of k-shapes and Genocchi numbers
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

def is_reducible2(p, k):
    # A method for determining k-reducibility of a k-shape
    rect_dim_list = k_rectangle_dimension_list(k) + k_rectangle_dimension_list(k-1)
    for (a, b) in rect_dim_list:
        if is_k_reducible_by_rectangle(p, k, (a,b)):
            return True
    return False

def is_reducible(ptn, k, method=2):
    """ A `k`-shape `ptn` is called *reducible* if there exists a `k`- or `k-1`-rectangle corresponding to both the `k`-row-shape and `k`-column-shape of `ptn`.  For a more rigorous definition, see Definition 3.7 of [genocchi]_.

    Note that this is different than the definition of a reducible partition!

    Given a `k`-shape `ptn` and a natural number `k`, returns True iff `ptn` is reducible.

    (Also, a `k`-shape is reducible iff it is not irreducible.)

    EXAMPLES::

        # The partition [3, 2, 1] has 3-row-shape [2, 2, 1] and 3-column-shape [2, 2, 1].  It is 3-reducible because there exists a 2x2-rectangle R in the 3-row-shape and the cells that make up R when viewed in the 3-column-shape form a 2x2-rectangle (you can't see it, but the 2's are switched here).
        sage: is_reducible(Partition([3, 2, 1]), k=3)
        True

        # In this example, no good rectangle can be found
        sage: is_reducible(Partition([5, 3, 2, 1, 1]), k=4)
        False
    """
    if method == 2:
        return is_reducible2(ptn, k)
    else:
        raise ValueError('Unknown reducibility method.')

def is_irreducible(s, k, method=2):
    """ A `k`-shape `ptn` is called *irreducible* if there does *not* exist a `k`- or `k-1`-rectangle corresponding to both the `k`-row-shape and `k`-column-shape of `ptn`.  For a more rigorous definition, see Definition 3.7 of [genocchi]_.

    Given a `k`-shape `ptn` and a natural number `k`, returns True iff `ptn` is irreducible.

    (Also, a `k`-shape is irreducible iff it is not reducible.)

    EXAMPLES::

        # The partition [3, 2, 1] has 3-row-shape [2, 2, 1] and 3-column-shape [2, 2, 1].  It is not 3-irreducible because there exists a 2x2-rectangle R in the 3-row-shape and the cells that make up R when viewed in the 3-column-shape form a 2x2-rectangle (you can't see it, but the 2's are switched here).
        sage: is_irreducible(Partition([3, 2, 1]), k=3)
        False

        # In this example, no good rectangle can be found, making it irreducible.
        sage: is_irreducible(Partition([5, 3, 2, 1, 1]), k=4)
        True
    """
    return not is_reducible(s, k, method)


############# GETTER FUNCS ############
def k_to_irreducible_k_shapes(k, method=2):
    """ Given a natural number `k`, return a list of all irreducible `k`-shapes.

    Note that the algorithm runs very slowly after `k=4` :(.

    EXAMPLES::

        sage: k_to_irreducible_k_shapes(3)
        [[], [1], [2, 1]]

    """
    bound = (k-1)*k/2
    n_bound = bound**2
    ptns = []
    for n in range(0, n_bound+1):
        ptns += Partitions(n, max_length=bound, max_part=bound)
        k_irr_k_shapes = [p for p in ptns if is_k_shape(p, k) and is_irreducible(p, k, method)]
    return k_irr_k_shapes
