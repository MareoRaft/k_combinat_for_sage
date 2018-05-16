#!/usr/bin/env sage
# A place to test my functions

from sage.all import *

from testing import *
from main import *
print('Sage loaded.  Testing...')


# test_right
sp = SkewPartition([[1], []])
c = right(sp, 0)
a(c, 0)

sp = SkewPartition([[2], []])
c = right(sp, 0)
a(c, 1)

sp = SkewPartition([[1, 1], []])
c = right(sp, 0)
a(c, 0)

sp = SkewPartition([[2], [1]])
c = right(sp, 0)
a(c, 1)

sp = SkewPartition([[2], [2]])
c = right(sp, 0)
a(c, None)

sp = SkewPartition([[5, 5, 4, 2, 2],  [4, 3, 2]])
c = right(sp, 2)
a(c, 3)


# test_left
sp = SkewPartition([[1], []])
c = left(sp, 0)
a(c, 0)

sp = SkewPartition([[2], []])
c = left(sp, 0)
a(c, 0)

sp = SkewPartition([[2], [1]])
c = left(sp, 0)
a(c, 1)

sp = SkewPartition([[2], [2]])
c = left(sp, 0)
a(c, None)

sp = SkewPartition([[5, 5, 4, 2, 2],  [4, 3, 2]])
c = left(sp, 2)
a(c, 2)


# test_top
sp = SkewPartition([[1], []])
c = top(sp, 0)
a(c, 0)

sp = SkewPartition([[2], []])
c = top(sp, 0)
a(c, 0)

sp = SkewPartition([[1, 1], []])
c = top(sp, 0)
a(c, 1)

sp = SkewPartition([[1, 1], [1]])
c = top(sp, 0)
a(c, 1)

sp = SkewPartition([[2], [2]])
c = top(sp, 0)
a(c, None)

sp = SkewPartition([[5, 5, 4, 2, 2],  [4, 3, 2]])
c = top(sp, 2)
a(c, 2)

sp = SkewPartition([[5, 5, 4, 2, 2],  [4, 3, 2]])
c = top(sp, 3)
a(c, 2)


# test_bottom
sp = SkewPartition([[1], []])
c = bottom(sp, 0)
a(c, 0)

sp = SkewPartition([[2], []])
c = bottom(sp, 0)
a(c, 0)

sp = SkewPartition([[1, 1], []])
c = bottom(sp, 0)
a(c, 0)

sp = SkewPartition([[1, 1], [1]])
c = bottom(sp, 0)
a(c, 1)

sp = SkewPartition([[2], [2]])
c = bottom(sp, 0)
a(c, None)

sp = SkewPartition([[5, 5, 4, 2, 2],  [4, 3, 2]])
c = bottom(sp, 2)
a(c, 2)

sp = SkewPartition([[5, 5, 4, 2, 2],  [4, 3, 2]])
c = bottom(sp, 3)
a(c, 1)


# test_bump_path_piece
sp = SkewPartition([[1], []])
top_row_index, is_end = bump_path_piece(sp, 0)
a(is_end, True)

sp = SkewPartition([[1, 1], []])
top_row_index, is_end = bump_path_piece(sp, 0)
a(is_end, True)

sp = SkewPartition([[2], []])
top_row_index, is_end = bump_path_piece(sp, 0)
a(is_end, True)

sp = SkewPartition([[3, 2, 1], [1]])
top_row_index, is_end = bump_path_piece(sp, 0)
a(is_end, False)
a(top_row_index, 2)

sp = SkewPartition([[6, 5, 3, 2, 2, 1], [2, 2]])
top_row_index, is_end = bump_path_piece(sp, 0)
a(is_end, False)
a(top_row_index, 3)

sp = SkewPartition([[6, 5, 3, 2, 2, 1], [2, 2]])
top_row_index, is_end = bump_path_piece(sp, 1, {0, 3})
a(top_row_index, 4)

sp = SkewPartition([[6, 5, 3, 2, 2, 1], [2, 2]])
top_row_index, is_end = bump_path_piece(sp, 2, {0, 1, 3, 4})
a(is_end, True)

sp = SkewPartition([[6, 5, 3, 2, 2, 1], [2, 2]])
top_row_index, is_end = bump_path_piece(sp, 5, {0, 1, 2, 3, 4})
a(is_end, True)


# test_bump_path
sp = SkewPartition([[1], []])
newly_blocked_rows = bump_path(sp, 0)
a(newly_blocked_rows, {0})

sp = SkewPartition([[1, 1], []])
newly_blocked_rows = bump_path(sp, 0)
a(newly_blocked_rows, {0})

sp = SkewPartition([[3, 2, 1], [1]])
newly_blocked_rows = bump_path(sp, 0)
a(newly_blocked_rows, {0, 2})

sp = SkewPartition([[3, 2, 1], [1]])
newly_blocked_rows = bump_path(sp, 1, {0, 2})
a(newly_blocked_rows, {1})

sp = SkewPartition([[6, 5, 3, 2, 2, 1], [2, 2]])
newly_blocked_rows = bump_path(sp, 0)
a(newly_blocked_rows, {0, 3})

sp = SkewPartition([[6, 5, 3, 2, 2, 1], [2, 2]])
newly_blocked_rows = bump_path(sp, 1, {0, 3})
a(newly_blocked_rows, {1, 4})

sp = SkewPartition([[6, 5, 3, 2, 2, 1], [2, 2]])
newly_blocked_rows = bump_path(sp, 2, {0, 3, 1, 4})
a(newly_blocked_rows, {2})


# test_skew_partition_to_selected_rows():
sp = SkewPartition([[], []])
selected_rows = skew_partition_to_selected_rows(sp)
a(selected_rows, [])

sp = SkewPartition([[1], []])
selected_rows = skew_partition_to_selected_rows(sp)
a(selected_rows, [0])

sp = SkewPartition([[3, 2, 1], [1]])
selected_rows = skew_partition_to_selected_rows(sp)
a(selected_rows, [0, 1])

sp = SkewPartition([[6, 5, 3, 2, 2, 1], [2, 2]])
selected_rows = skew_partition_to_selected_rows(sp)
a(selected_rows, [0, 1, 2, 5])


# test_selected_rows_to_root_ideal
root_ideal = selected_rows_to_root_ideal(1, [0])
a(root_ideal, [])

root_ideal = selected_rows_to_root_ideal(2, [0, 1])
a(root_ideal, [])

root_ideal = selected_rows_to_root_ideal(2, [0])
a(root_ideal, [(0,1)])

root_ideal = selected_rows_to_root_ideal(3, [0, 2])
a(root_ideal, [(0,1), (0,2)])

root_ideal = selected_rows_to_root_ideal(4, [0, 2])
a(root_ideal, [(0,1), (0,2), (0,3), (1,3)])

root_ideal = selected_rows_to_root_ideal(5, [0, 1, 4])
a(root_ideal, [(0,2), (0,3), (0,4), (1,3), (1,4)])

root_ideal = selected_rows_to_root_ideal(5, [0, 1])
a(root_ideal, [(0,2), (0,3), (0,4), (1,3), (1,4), (2,4)])


# test_skew_partition_to_root_ideal
sp = SkewPartition([[6, 5, 3, 2, 2, 1], [2, 2]])
ri = skew_partition_to_root_ideal(sp)
a(ri, [(0,3), (0,4), (0,5), (1,4), (1,5)])


# test_k_row_lengths
ptn = Partition([4, 4, 4, 3, 2])
rs = k_row_lengths(ptn, 2)
a(rs, [0, 1, 1, 1, 2])


# test_has_rectangle


# test_has_k_rectangle


#


# TESTS:
#     # empty skew
#     sage: sp = SkewPartition2([[], []])
#     sage: sp.is_skew_linked_diagram()
#     True
#     # # valid row shape but invalid col shape
#     # sage: sp = SkewPartition2([[3, 2], [1, 0]])
#     # sage: sp.is_skew_linked_diagram()
#     # False
#     # sage: aFalse


# ALL DONE!
print('Testing completed successfully!')
