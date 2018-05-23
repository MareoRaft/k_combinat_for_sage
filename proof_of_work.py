#!/usr/bin/env sage
# -*- coding: utf-8 -*-
""" This file demonstrates that I have completed each requested task (or the functionality already existed in sage) by showing an example for each one. """
from sage.all import *
from main import *
print('Sage loaded.  Executing proof of work...')


# 1. Given a skew-shape, output the pair (row-shape, col-shape).
sp = SkewPartition([[7, 7, 3, 3, 2, 2],  [6, 4, 1, 1]])
pair = (sp.row_lengths(), sp.column_lengths())
assert pair == ([1, 3, 2, 2, 2, 2], [2, 4, 2, 0, 1, 1, 2])
# note that there is no "pair" method like sp.row_column_lengths_pair.


# 2. Given a (row-shape, col-shape) pair, output the skew-shape or error.
sp = row_col_to_skew_partition([4, 2, 2, 1, 1, 1, 1], [3, 2, 2, 1, 1, 1, 1, 1])
assert sp.to_list() == [[8, 4, 3, 2, 1, 1, 1], [4, 2, 1, 1]]
# Actually, i think this already exists as SkewPartitions().from_row_and_column_length, but I can't get that to run


# 3. Given a skew-linked diagram (or maybed just skew-shape) and k, detect if it is a k-boundary.
sp = SkewPartition([[5, 5, 4, 2, 2],  [4, 3, 2]])
k = 3
assert is_k_boundary(sp, k) == True


# 4. Given k, output a list of all k-irreducible partitions
irr_ptns = get_k_irreducible_partitions(3)
assert irr_ptns == [Partition([]), Partition([1]), Partition([1, 1]), Partition([2]), Partition([2, 1]), Partition([2, 1, 1])]
# or for efficiency, get lists instead of actual Partition objects
irr_ptns = get_k_irreducible_partition_lists(3)
assert irr_ptns == [[], [1], [1, 1], [2], [2, 1], [2, 1, 1]]


# 5. Given skew-linked diagram, generate root ideal (in progress)
sp = SkewPartition([[6, 5, 3, 2, 2, 1], [2, 2]])
root_ideal = skew_partition_to_root_ideal(sp)
assert root_ideal == [(0,3), (0,4), (0,5), (1,4), (1,5)]
# note that all indecis are 0-based


# 6. Given a partition λ, detect if λ = λ'
p = Partition([6, 4, 4, 3, 1, 1])
assert is_symmetric(p)


# 7. Get the first 20 terms of the sequence (a_n) where a_n is the number of pairs (λ, λ) where λ has size n and there exists a skew-partition whose row and col shapes are λ
seq = sequence(n_to_number_of_linked_partition_self_pairs)
assert seq == [1, 1, 1, 2, 3, 4, 4, 7, 9, 13, 12, 20, 24, 32, 31, 50, 55, 74, 76, 109]
# To get only the first 10 terms, use sequence(n_to_number_of_linked_partition_self_pairs, num_terms=10).
# To print out a sequence, use print_sequence(n_to_number_of_linked_partition_self_pairs).


# 8. Find some more interesting sequences:
# Get the sequence (a_n) where where a_n is the number of 1-shapes of size n.
seq = sequence(lambda n: n_to_num_k_shapes(n, k=1))
assert seq == [1, 1, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0] # A010054
seq = sequence(lambda n: n_to_num_k_shapes(n, k=2))
assert seq == [1, 1, 2, 1, 2, 1, 3, 2, 1, 2, 3, 2, 3, 2, 2, 1, 5, 3, 2, 2]
seq = sequence(lambda n: n_to_num_k_shapes(n, k=3))
assert seq == [1, 1, 2, 3, 3, 3, 5, 5, 5, 8, 6, 6, 10, 9, 11, 10, 9, 13, 15, 13]
# n to number of k-shapes of size n (for any k between 1 and n-1)
seq = sequence(n_to_num_k_shapes)
assert seq == [0, 0, 0, 1, 3, 5, 9, 13, 20, 28, 40, 54, 75, 99, 133, 174, 229, 295, 383, 488]
# n to number of self-conjugate k-skews
seq = sequence(lambda n: n_to_num_self_conjugate_k_skews(n, k=0))
assert seq == [1, 1, 0, 1, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 4, 5, 5, 5, 6] # A000700
seq = sequence(lambda n: n_to_num_self_conjugate_k_skews(n, 2))
assert seq == [1, 1, 0, 1, 0, 1, 1, 0, 1, 0, 1, 0, 1, 0, 0, 1, 1, 1, 0, 0]
# TODO: number of k-skews of size n.


# 9. Given k, find all k-irreducible-k-shapes.
ptns = get_k_irreducible_k_shapes(2)
assert ptns == [[], [1], [2, 1], [3, 2, 1]]



# ALL DONE!
print('Proof of work completed successfully!')
