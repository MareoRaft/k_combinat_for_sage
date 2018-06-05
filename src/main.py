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

def n_to_self_conjugate_k_shape_boundaries(n, k):
    k_shape_boundaries = n_to_k_shape_boundaries(n, k)
    return [ks for ks in k_shape_boundaries if ks == ks.conjugate()]

def n_to_num_self_conjugate_k_shape_boundaries(n, k):
    return len(n_to_self_conjugate_k_shape_boundaries(n, k))



