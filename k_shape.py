# -*- coding: utf-8 -*-
from sage.all import *
from partition import *
import skew_partition as SP

# HELPERS
def k_rectangle_dimension_list(k):
    return [(k-i+1, i) for i in range(1, k+1)]

# kShape verifier
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
        return SP.is_linked(k_bdy)

#kShape stuff
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

def is_k_reducible_by_rectangle(p, k, (a,b)):
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

def is_k_reducible2(p, k):
    rect_dim_list = k_rectangle_dimension_list(k) + k_rectangle_dimension_list(k-1)
    for (a, b) in rect_dim_list:
        if is_k_reducible_by_rectangle(p, k, (a,b)):
            return True
    return False

def is_k_reducible(s, k, method=2):
    if method == 2:
        return is_k_reducible2(s, k)
    else:
        raise ValueError('Unknown reducibility method.')

def is_k_irreducible(s, k, method=2):
    """ A k-shape is called __k-irreducible__ if it is not k-reducible. """
    return not is_k_reducible(s, k, method)


############# GETTER FUNCS ############
def get_k_irreducible_k_shapes(k, method=2):
    # The k-row-shape has at most k rows of length 0, k-1 rows of length 1, ..., 0 rows of length k.  And 0 rows of length greater than k.  Hence the k-row-shape has an upper bound of k*(k-1)/2 rows.  The same goes for the k-col-shape.
    bound = (k-1)*k/2
    n_bound = bound**2
    ptns = []
    for n in range(0, n_bound+1):
        ptns += Partitions(n, max_length=bound, max_part=bound)
        k_irr_k_shapes = [p for p in ptns if is_k_shape(p, k) and is_k_irreducible(p, k, method)]
    return k_irr_k_shapes
