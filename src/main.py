#!/usr/bin/env sage
# -*- coding: utf-8 -*-
from sage.all import *
from partition import *
import partition as P
from skew_partition import *
import skew_partition as SP
from k_shape import *
import k_shape as kS
from root_ideal import *
import root_ideal as RI



# MAIN
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

def n_to_k_shapes(n, k):
    """ Given n, find all partitions of size n that are k-shapes. """
    return [ptn for ptn in Partitions(n) if is_k_shape(ptn, k)]

def n_to_num_k_shapes(n, k):
    return len(n_to_k_shapes(n, k))

def n_to_k_shape_boundaries(n, k):
    """ Given n, find all k-boundaries of all k-shapes of size n. """
    return [ptn.k_boundary(k) for ptn in Partitions(n) if is_k_shape(ptn, k)]

def n_to_symmetric_k_shape_boundaries(n, k):
    k_shape_boundaries = n_to_k_shape_boundaries(n, k)
    return [ks for ks in k_shape_boundaries if kS.is_symmetric(ks)]

def n_to_num_symmetric_k_shape_boundaries(n, k):
    return len(n_to_symmetric_k_shape_boundaries(n, k))






## Compositional hall-littlewood polynomials and more
# def hall_littlewood_vertex_function_base(m, input_, base_ring=QQ):
#     assert m == ZZ(m)
#     # add 't' to ring and bless t variable
#     base_ring = base_ring['t']
#     t = base_ring.gen()
#     sym = SymmetricFunctions(base_ring)
#     h = sym.h()
#     e = sym.e()
#     def summand(i, j, input_):
#         h_a = h([m + i + j])
#         step1 = input_.skew_by(h[j])
#         step2 = step1.skew_by(e[i])
#         step3 = step2 * (-1)**i * t**j * h_a
#         return step3
#     i_max = 99
#     j_max = 99
#     return sum(summand(i, j, input_) for i in range(0, i_max+1) for j in range(0, j_max+1))
# def hall_littlewood_vertex_function(ptn, input_, base_ring=QQ):
#     # base cases
#     if isinstance(ptn, int):
#         return hall_littlewood_vertex_function_base(ptn, input_, base_ring)
#     if len(ptn) == 0:
#         return input_
#     if len(ptn) == 1:
#         return hall_littlewood_vertex_function_base(ptn[0], input_, base_ring)
#     # inductive step
#     for part in reversed(ptn):
#         input_ = hall_littlewood_vertex_function_base(part, input_, base_ring)
#     return input_
# # def hall_littlewood_vertex_operator(ptn, base_ring=QQ):
# #     """
# #     base_ring: the base ring to build the SymmetricFunctions upon.

# #     see p.14 of [cat]_.
# #     """
# #     return lambda input_: hall_littlewood_vertex_function(ptn, input_, base_ring)
# class HallLittlewoodVertexOperator:
#     """
#     base_ring: the base ring to build the SymmetricFunctions upon.

#     see p.14 of [cat]_.
#     """
#     def __init__(self, ptn, base_ring=QQ):
#         self.ptn = ptn
#         self.base_ring = base_ring

#     def __call__(self, input_):
#         # hall_littlewood_vertex_operator
#         return hall_littlewood_vertex_function(self.ptn, input_, self.base_ring)









