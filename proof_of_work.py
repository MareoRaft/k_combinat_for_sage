#!/usr/bin/env sage
# -*- coding: utf-8 -*-
""" This file demonstrates that I have completed each requested task (or the functionality already existed in sage) by showing an example for each one. """
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


# 3. Given a skew-linked diagram (or maybed just skew-shape) and k, detect if it is a k-boundary.
sp = SkewPartition([[5, 5, 4, 2, 2],  [4, 3, 2]])
k = 3
assert is_k_boundary(sp, k) == True


# 4. Given k, output a list of all k-irreducible partitions
irr_partitions = get_k_irreducible_partition_lists(3)
assert irr_partitions == [[], [1], [1, 1], [2], [2, 1], [2, 1, 1]]


# 5. Given skew-linked diagram, generate root ideal (in progress)
# sp = SkewPartition([[6, 5, 3, 2, 2, 1] [2, 2]])
# selected_rows = skew_partition_to_selected_rows(sp)
# root_ideal = selected_rows_to_root_ideal(selected_rows)
# assert root_ideal == [(1,4), (1,5), (1,6), (2,5), (2,6)]


# 6. Given a partition λ, detect if λ = λ'
p = Partition([6, 4, 4, 3, 1, 1])
assert is_symmetric(p)




# ALL DONE!
print('Proof of work completed successfully!')
