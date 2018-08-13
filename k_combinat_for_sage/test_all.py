#!/usr/bin/env sage
# use `sage --python -m pdb test_all.py` for the debugger
# A place to test my functions
from __future__ import print_function
import time

from sage.all import *
print('Sage loaded.  Now loading local modules...')
from testing import *
from all import *
from all import _is_sequence
from strong_marked_tableau import __go_to_ribbon_head
start_time = time.time()
print('Modules loaded.  Testing...')



# test summands / terms
s = SymmetricFunctions(QQ).s()
a(set((s[1] + 2 + 3*s[34]).terms()), set([2*s[[]], s[1], 3*s[34]]))


# test k size
a(k_size([2, 1, 1], 1), 2)
a(k_size([2, 1, 1], 2), 3)
a(k_size([2, 1, 1], 3), 3)
a(k_size([2, 1, 1], 4), 4)


# test is sequence
a(_is_sequence([1, 3, 2]), True)


# test_SP_is_symmetric
sp = SkewPartition([[], []])
b = skew_partition.is_symmetric(sp)
a(b, True)

sp = SkewPartition([[1], []])
b = skew_partition.is_symmetric(sp)
a(b, True)

sp = SkewPartition([[4,3,3,1], [1]])
b = skew_partition.is_symmetric(sp)
a(b, True)

sp = SkewPartition([[4,3,3,1], [1,1]])
b = skew_partition.is_symmetric(sp)
a(b, False)

sp = SkewPartition([[5,3,3,1], [2,2]])
b = skew_partition.is_symmetric(sp)
a(b, False)



# test_right
sp = SkewPartition([[1], []])
c = skew_partition.right(sp, 0)
a(c, 0)

sp = SkewPartition([[2], []])
c = skew_partition.right(sp, 0)
a(c, 1)

sp = SkewPartition([[1, 1], []])
c = skew_partition.right(sp, 0)
a(c, 0)

sp = SkewPartition([[2], [1]])
c = skew_partition.right(sp, 0)
a(c, 1)

sp = SkewPartition([[2], [2]])
c = skew_partition.right(sp, 0)
a(c, None)

sp = SkewPartition([[5, 5, 4, 2, 2],  [4, 3, 2]])
c = skew_partition.right(sp, 2)
a(c, 3)


# test_left
sp = SkewPartition([[1], []])
c = skew_partition.left(sp, 0)
a(c, 0)

sp = SkewPartition([[2], []])
c = skew_partition.left(sp, 0)
a(c, 0)

sp = SkewPartition([[2], [1]])
c = skew_partition.left(sp, 0)
a(c, 1)

sp = SkewPartition([[2], [2]])
c = skew_partition.left(sp, 0)
a(c, None)

sp = SkewPartition([[5, 5, 4, 2, 2],  [4, 3, 2]])
c = skew_partition.left(sp, 2)
a(c, 2)


# test_top
sp = SkewPartition([[1], []])
c = skew_partition.top(sp, 0)
a(c, 0)

sp = SkewPartition([[2], []])
c = skew_partition.top(sp, 0)
a(c, 0)

sp = SkewPartition([[1, 1], []])
c = skew_partition.top(sp, 0)
a(c, 1)

sp = SkewPartition([[1, 1], [1]])
c = skew_partition.top(sp, 0)
a(c, 1)

sp = SkewPartition([[2], [2]])
c = skew_partition.top(sp, 0)
a(c, None)

sp = SkewPartition([[5, 5, 4, 2, 2],  [4, 3, 2]])
c = skew_partition.top(sp, 2)
a(c, 2)

sp = SkewPartition([[5, 5, 4, 2, 2],  [4, 3, 2]])
c = skew_partition.top(sp, 3)
a(c, 2)


# test_bottom
sp = SkewPartition([[1], []])
c = skew_partition.bottom(sp, 0)
a(c, 0)

sp = SkewPartition([[2], []])
c = skew_partition.bottom(sp, 0)
a(c, 0)

sp = SkewPartition([[1, 1], []])
c = skew_partition.bottom(sp, 0)
a(c, 0)

sp = SkewPartition([[1, 1], [1]])
c = skew_partition.bottom(sp, 0)
a(c, 1)

sp = SkewPartition([[2], [2]])
c = skew_partition.bottom(sp, 0)
a(c, None)

sp = SkewPartition([[5, 5, 4, 2, 2],  [4, 3, 2]])
c = skew_partition.bottom(sp, 2)
a(c, 2)

sp = SkewPartition([[5, 5, 4, 2, 2],  [4, 3, 2]])
c = skew_partition.bottom(sp, 3)
a(c, 1)


# test_bounce_path_piece
# sp = SkewPartition([[1], []])
# top_row_index, is_end = bounce_path_piece(sp, 0)
# a(is_end, True)

# sp = SkewPartition([[1, 1], []])
# top_row_index, is_end = bounce_path_piece(sp, 0)
# a(is_end, True)

# sp = SkewPartition([[2], []])
# top_row_index, is_end = bounce_path_piece(sp, 0)
# a(is_end, True)

# sp = SkewPartition([[3, 2, 1], [1]])
# top_row_index, is_end = bounce_path_piece(sp, 0)
# a(is_end, False)
# a(top_row_index, 2)

# sp = SkewPartition([[6, 5, 3, 2, 2, 1], [2, 2]])
# top_row_index, is_end = bounce_path_piece(sp, 0)
# a(is_end, False)
# a(top_row_index, 3)

# sp = SkewPartition([[6, 5, 3, 2, 2, 1], [2, 2]])
# top_row_index, is_end = bounce_path_piece(sp, 1, {0, 3})
# a(top_row_index, 4)

# sp = SkewPartition([[6, 5, 3, 2, 2, 1], [2, 2]])
# top_row_index, is_end = bounce_path_piece(sp, 2, {0, 1, 3, 4})
# a(is_end, True)

# sp = SkewPartition([[6, 5, 3, 2, 2, 1], [2, 2]])
# top_row_index, is_end = bounce_path_piece(sp, 5, {0, 1, 2, 3, 4})
# a(is_end, True)


# test_bounce_path
# sp = SkewPartition([[1], []])
# newly_blocked_rows = bounce_path(sp, 0)
# a(newly_blocked_rows, {0})

# sp = SkewPartition([[1, 1], []])
# newly_blocked_rows = bounce_path(sp, 0)
# a(newly_blocked_rows, {0})

# sp = SkewPartition([[3, 2, 1], [1]])
# newly_blocked_rows = bounce_path(sp, 0)
# a(newly_blocked_rows, {0, 2})

# sp = SkewPartition([[3, 2, 1], [1]])
# newly_blocked_rows = bounce_path(sp, 1, {0, 2})
# a(newly_blocked_rows, {1})

# sp = SkewPartition([[6, 5, 3, 2, 2, 1], [2, 2]])
# newly_blocked_rows = bounce_path(sp, 0)
# a(newly_blocked_rows, {0, 3})

# sp = SkewPartition([[6, 5, 3, 2, 2, 1], [2, 2]])
# newly_blocked_rows = bounce_path(sp, 1, {0, 3})
# a(newly_blocked_rows, {1, 4})

# sp = SkewPartition([[6, 5, 3, 2, 2, 1], [2, 2]])
# newly_blocked_rows = bounce_path(sp, 2, {0, 3, 1, 4})
# a(newly_blocked_rows, {2})


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


# test_skew_partition_to_removable_roots
sp = SkewPartition([[], []])
a(skew_partition_to_removable_roots(sp, type='max'), [])

sp = SkewPartition([[3, 2, 1], [1]])
a(skew_partition_to_removable_roots(sp, type='max'), [(0,2)])

sp = SkewPartition([[6, 5, 3, 2, 2, 1], [2, 2]])
a(skew_partition_to_removable_roots(sp, type='max'), [(0,3), (1,4)])



# test_selected_rows_to_maximum_root_ideal
ri = selected_rows_to_maximum_root_ideal(1, [0])
a(ri, [])

ri = selected_rows_to_maximum_root_ideal(2, [0, 1])
a(ri, [])

ri = selected_rows_to_maximum_root_ideal(2, [0])
a(ri, [(0,1)])

ri = selected_rows_to_maximum_root_ideal(3, [0, 2])
a(ri, [(0,1), (0,2)])

ri = selected_rows_to_maximum_root_ideal(4, [0, 2])
a(ri, [(0,1), (0,2), (0,3), (1,3)])

ri = selected_rows_to_maximum_root_ideal(5, [0, 1, 4])
a(ri, [(0,2), (0,3), (0,4), (1,3), (1,4)])

ri = selected_rows_to_maximum_root_ideal(5, [0, 1])
a(ri, [(0,2), (0,3), (0,4), (1,3), (1,4), (2,4)])


# test_removable_roots_to_partition
rr = []
n = 3
a(removable_roots_to_partition(rr, n), [])

rr = [(0, 0)]
n = 1
a(removable_roots_to_partition(rr, n), [1])

rr = [(0, 0)]
n = 2
a(removable_roots_to_partition(rr, n), [2])

rr = [(1, 1)]
n = 2
a(removable_roots_to_partition(rr, n), [1, 1])

rr = [(1, 1)]
n = 2
a(removable_roots_to_partition(rr, n), [1, 1])

rr = [(0, 3), (1, 4)]
n = 6
a(removable_roots_to_partition(rr, n), [3, 2])


# test_skew_RIS.init_from_partition
RIS = RootIdeals()
sp = SkewPartition([[6, 5, 3, 2, 2, 1], [2, 2]])
ri = RIS.init_from_skew_partition(sp, type='max')
a(ri, [(0,3), (0,4), (0,5), (1,4), (1,5)])

sp = SkewPartition([[4, 3, 2, 2, 1, 1], [3, 2, 1, 1]])
ri = RIS.init_from_skew_partition(sp, type='max')
a(ri, [(0,1), (0,2), (0,3), (0,4), (0,5), (1,2), (1,3), (1,4), (1,5), (2,4), (2,5), (3,5)])

sp = SkewPartition([[4, 3, 2, 2, 1, 1], [3, 2, 1, 1]])
ri = RIS.init_from_skew_partition(sp, type='min')
a(ri, [(0,1), (0,2), (0,3), (0,4), (0,5), (1,3), (1,4), (1,5), (2,4), (2,5), (3,5)])


# test_removable_roots_to_root_ideal
rr = [(0,1), (2,2)]
n = 4
a(RIS.init_from_removable_roots(rr, n), [(0,1), (0,2), (0,3), (1,2), (1,3), (2,2), (2,3)])


# test_down
# Rightmost example 2.4 of SKEW-LINKED CATALAN FUNCTIONS AND k-SCHUR POSITIVITY
ri = RootIdeal([(0,1), (0,2), (0,3), (0,4), (0,5), (1,4), (1,5), (2,4), (2,5), (3,4), (3,5)])
a(ri.down(0), 1)
a(ri.down(1), 4)
a(ri.down(2), 4)
a(ri.down(3), 4)
a(ri.down(4), None)
a(ri.down(5), None)


# test_down_path
# Rightmost example 2.4 of SKEW-LINKED CATALAN FUNCTIONS AND k-SCHUR POSITIVITY
ri = RootIdeal([(0,1), (0,2), (0,3), (0,4), (0,5), (1,4), (1,5), (2,4), (2,5), (3,4), (3,5)])
a(ri.down_path(0), [0, 1, 4])
a(ri.down_path(1), [1, 4])
a(ri.down_path(2), [2, 4])
a(ri.down_path(3), [3, 4])
a(ri.down_path(4), [4])
a(ri.down_path(5), [5])


# test_down_path_column_lengths_part
# Rightmost example 2.4 of SKEW-LINKED CATALAN FUNCTIONS AND k-SCHUR POSITIVITY
ri = RootIdeal([(0,1), (0,2), (0,3), (0,4), (0,5), (1,4), (1,5), (2,4), (2,5), (3,4), (3,5)])
ptn = [7, 6, 5, 2, 2, 2]
a(ri.down_path_column_lengths_part(ptn, 0), 15)
a(ri.down_path_column_lengths_part(ptn, 1), 8)
a(ri.down_path_column_lengths_part(ptn, 2), 7)
a(ri.down_path_column_lengths_part(ptn, 3), 4)
a(ri.down_path_column_lengths_part(ptn, 4), 2)
a(ri.down_path_column_lengths_part(ptn, 5), 2)


# test_down_path_column_lengths
# Rightmost example 2.4 of SKEW-LINKED CATALAN FUNCTIONS AND k-SCHUR POSITIVITY
ri = RootIdeal([(0,1), (0,2), (0,3), (0,4), (0,5), (1,4), (1,5), (2,4), (2,5), (3,4), (3,5)])
ptn = [7, 6, 5, 2, 2, 2]
a(ri.down_path_column_lengths(ptn), [15, 7, 4, 2])


# test_root_ideal_to_partition
ri = RootIdeal([])
a(ri.to_partition(), [])

# Rightmost example 2.4 of SKEW-LINKED CATALAN FUNCTIONS AND k-SCHUR POSITIVITY
ri = RootIdeal([(0,1), (0,2), (0,3), (0,4), (0,5), (1,4), (1,5), (2,4), (2,5), (3,4), (3,5)])
a(ri.to_partition(), [5, 2, 2, 2])


# test_RIS.init_from_partition
p = Partition([])
n = 5
a(RIS.init_from_partition(p, n), [])

# Rightmost example 2.4 of SKEW-LINKED CATALAN FUNCTIONS AND k-SCHUR POSITIVITY
p = Partition([5, 2, 2, 2])
n = 6
a(RIS.init_from_partition(p, n), [(0,1), (0,2), (0,3), (0,4), (0,5), (1,4), (1,5), (2,4), (2,5), (3,4), (3,5)])


# test_is_strict
ri = RootIdeal([])
a(ri.is_strict(), True)

# Rightmost example 2.4 of SKEW-LINKED CATALAN FUNCTIONS AND k-SCHUR POSITIVITY
ri = RootIdeal([(0,1), (0,2), (0,3), (0,4), (0,5), (1,4), (1,5), (2,4), (2,5), (3,4), (3,5)])
a(ri.is_strict(), False)

ri = RootIdeal([(0,1), (0,2), (0,3), (0,4), (0,5)])
a(ri.is_strict(), True)

ri = RootIdeal([(0,1), (0,2), (0,3), (0,4), (0,5), (1,2), (1,3), (1,4), (1,5), (2,3), (2,4), (2,5), (3,4), (3,5)])
a(ri.is_strict(), True)

ri = RootIdeal([(0,1), (0,2), (0,3), (0,4), (0,5), (1,2), (1,3), (1,4), (1,5), (2,3), (2,4), (2,5), (3,4), (3,5), (4,5)])
a(ri.is_strict(), True)


# test_boundary
p = Partition([1])
a(boundary(p), [(1,0), (1,1), (0,1)])

p = Partition([2, 1])
a(boundary(p), [(2,0), (2,1), (1,1), (1,2), (0,2)])

p = Partition([3, 1])
a(boundary(p), [(3,0), (3,1), (2,1), (1,1), (1,2), (0,2)])

p = Partition([2, 1, 1])
a(boundary(p), [(2,0), (2,1), (1,1), (1,2), (1,3), (0,3)])

p = Partition([7, 4, 3, 3, 2, 1, 1])
a(boundary(p), [(7,0), (7,1), (6,1), (5,1), (4,1), (4,2), (3,2), (3,3), (3,4), (2,4), (2,5), (1,5), (1,6), (1,7), (0,7)])


# test_k_rim
p = Partition([1])
k = 0
a(k_rim(p, k), [(1,0), (1,1), (0,1)])

p = Partition([3, 1])
k = 0
a(k_rim(p, k), [(3,0), (3,1), (2,1), (1,1), (1,2), (0,2)])

p = Partition([3, 1])
k = 1
a(k_rim(p, k), [(3,0), (2,0), (2,1), (1,1), (0,1), (0,2)])

p = Partition([3, 1])
k = 2
a(k_rim(p, k), [(3,0), (2,0), (1,0), (1,1), (0,1), (0,2)])

p = Partition([3, 1])
k = 3
a(k_rim(p, k), [(3,0), (2,0), (1,0), (1,1), (0,1), (0,2)])

p = Partition([3, 1])
k = 4
a(k_rim(p, k), [(3,0), (2,0), (1,0), (0,0), (0,1), (0,2)])

p = Partition([7, 4, 3, 3, 2, 1, 1])
k = 0
a(k_rim(p, k), [(7,0), (7,1), (6,1), (5,1), (4,1), (4,2), (3,2), (3,3), (3,4), (2,4), (2,5), (1,5), (1,6), (1,7), (0,7)])

p = Partition([7, 4, 3, 3, 2, 1, 1])
k = 1
a(k_rim(p, k), [(7,0), (6,0), (6,1), (5,1), (4,1), (3,1), (3,2), (3,3), (2,3), (2,4), (1,4), (1,5), (1,6), (0,6), (0,7)])




# test_k_row_lengths
ptn = Partition([4, 4, 4, 3, 2])
rs = k_row_lengths(ptn, 2)
a(rs, [0, 1, 1, 1, 2])


# test_has_rectangle
p = Partition([1, 1, 1])
h = 4
w = 1
a(has_rectangle(p, h, w), False)

p = Partition([1, 1, 1])
h = 3
w = 1
a(has_rectangle(p, h, w), True)

p = Partition([1, 1, 1])
h = 2
w = 1
a(has_rectangle(p, h, w), True)

p = Partition([1, 1, 1])
h = 1
w = 2
a(has_rectangle(p, h, w), False)

p = Partition([3])
h = 1
w = 3
a(has_rectangle(p, h, w), True)

p = Partition([3])
h = 1
w = 2
a(has_rectangle(p, h, w), False)

p = Partition([3])
h = 2
w = 3
a(has_rectangle(p, h, w), False)

p = Partition([5, 5, 4, 2, 2])
h = 2
w = 2
a(has_rectangle(p, h, w), True)

p = Partition([5, 5, 4, 2, 2])
h = 3
w = 2
a(has_rectangle(p, h, w), False)

p = Partition([5, 5, 4, 2, 2])
h = 2
w = 5
a(has_rectangle(p, h, w), True)


# test_has_k_rectangle
p = Partition([1])
k = 1
a(has_k_rectangle(p, k), True)

p = Partition([1])
k = 2
a(has_k_rectangle(p, k), False)

p = Partition([1, 1, 1])
k = 3
a(has_k_rectangle(p, k), True)

p = Partition([1, 1, 1])
k = 2
a(has_k_rectangle(p, k), True)

p = Partition([1, 1, 1])
k = 4
a(has_k_rectangle(p, k), False)

p = Partition([3])
k = 3
a(has_k_rectangle(p, k), True)

# p = Partition([3]) # DON"T KNOW
# k = 2
# a(has_k_rectangle(p, k), ?)

p = Partition([3])
k = 4
a(has_k_rectangle(p, k), False)

p = Partition([5, 5, 4, 2, 2])
k = 7
a(has_k_rectangle(p, k), False)

p = Partition([5, 5, 4, 2, 2])
k = 6
a(has_k_rectangle(p, k), True)

p = Partition([5, 5, 4, 2, 2])
k = 5
a(has_k_rectangle(p, k), True)

p = Partition([5, 5, 4, 2, 2])
k = 4
a(has_k_rectangle(p, k), True)

p = Partition([5, 5, 4, 2, 2])
k = 3
a(has_k_rectangle(p, k), True)

p = Partition([5, 5, 4, 2, 2])
k = 2
a(has_k_rectangle(p, k), True)


# test_partition.next
min_ = []
max_ = []
p = []
a(partition.next(p, min=min_, max=max_), func=is_False_or_None)

min_ = []
max_ = [1]
a(partition.next([], min=min_, max=max_), [1])
a(partition.next([1], min=min_, max=max_), func=is_False_or_None)

min_ = []
max_ = [1, 1]
a(partition.next([], min=min_, max=max_), [1])
a(partition.next([1], min=min_, max=max_), [1, 1])
a(partition.next([1, 1], min=min_, max=max_), func=is_False_or_None)

min_ = []
max_ = [2]
a(partition.next([], min=min_, max=max_), [1])
a(partition.next([1], min=min_, max=max_), [2])
a(partition.next([2], min=min_, max=max_), func=is_False_or_None)

min_ = []
max_ = [2, 1]
a(partition.next([], min=min_, max=max_), [1])
a(partition.next([1], min=min_, max=max_), [1, 1])
a(partition.next([1, 1], min=min_, max=max_), [2])
a(partition.next([2], min=min_, max=max_), [2, 1])
a(partition.next([2, 1], min=min_, max=max_), func=is_False_or_None)

min_ = [1, 1]
max_ = [3, 2, 1]
a(partition.next([1, 1], min=min_, max=max_), [1, 1, 1])
a(partition.next([1, 1, 1], min=min_, max=max_), [2, 1])


# test_RootIdeal_next
ri = RootIdeal([], n=4)
a(ri.next(), [(0,3)])

max_ri = RootIdeal([(0,0), (0,1), (1,0), (1,1)])
min_ri = RootIdeal([])
ri = RootIdeal([(0,1), (1,1)])
a(ri.next(min=min_ri, max=max_ri), [(0,0), (0,1)])

max_ri = RootIdeal([(0,1), (0,2), (0,3), (1,2), (1,3), (2,3)])
min_ri = RootIdeal([(0,2), (0,3)])
ri = RootIdeal([(0,1), (0,2), (0,3)])
a(ri.next(min=min_ri, max=max_ri), [(0,1), (0,2), (0,3), (1,3)])


# test_RIS.init_all_from_skew_partition
sp = SkewPartition([[4, 2, 1, 1], [2, 1]])
correct_ris = [
	RootIdeal([(0,1), (0,2), (0,3), (1,2), (1,3)]),
	RootIdeal([(0,1), (0,2), (0,3), (1,3)]),
]
a(set(RIS.init_all_from_skew_partition(sp)), set(correct_ris))

sp = SkewPartition([[4, 3, 2, 2, 1, 1], [3, 2, 1, 1]])
correct_ris = [
	RootIdeal([(0,1), (0,2), (0,3), (0,4), (0,5), (1,2), (1,3), (1,4), (1,5), (2,4), (2,5), (3,5)]),
	RootIdeal([(0,1), (0,2), (0,3), (0,4), (0,5), (1,3), (1,4), (1,5), (2,4), (2,5), (3,5)]),
]
a(set(RIS.init_all_from_skew_partition(sp)), set(correct_ris))


# test_get_k_rectangles
out = set(get_k_rectangles(0))
a(out, set(Partition([])))

out = set(get_k_rectangles(1))
a(out, set([Partition([1])]))

out = set(get_k_rectangles(2))
a(out, set([Partition([1, 1]), Partition([2])]))

out = set(get_k_rectangles(3))
a(out, set([Partition([1, 1, 1]), Partition([2, 2]), Partition([3])]))

out = set(get_k_rectangles(4))
a(out, set([Partition([1, 1, 1, 1]), Partition([2, 2, 2]), Partition([3, 3]), Partition([4])]))


# test_h_bounds
p = Partition([10, 7, 4, 2, 2, 2, 1, 1, 1, 1])
k = 4
i = 1
a(h_bounds(p, k, i), (3, 10))

p = Partition([10, 7, 4, 2, 2, 2, 1, 1, 1, 1])
k = 4
i = 2
a(h_bounds(p, k, i), (2, 3))

p = Partition([10, 7, 4, 2, 2, 2, 1, 1, 1, 1])
k = 4
i = 3
a(h_bounds(p, k, i), (0, 2))


# test_v_bounds
p = Partition([10, 7, 4, 2, 2, 2, 1, 1, 1, 1])
k = 4
i = 1
a(v_bounds(p, k, i), (2, 10))

p = Partition([10, 7, 4, 2, 2, 2, 1, 1, 1, 1])
k = 4
i = 2
a(v_bounds(p, k, i), (2, 2))

p = Partition([10, 7, 4, 2, 2, 2, 1, 1, 1, 1])
k = 4
i = 3
a(v_bounds(p, k, i), (1, 2))

p = Partition([10, 7, 4, 2, 2, 2, 1, 1, 1, 1])
k = 4
i = 4
a(v_bounds(p, k, i), (0, 1))


# test_k_shape.is_k_reducible_by_rectangle
s = Partition([1])
k = 1
(w,h) = (1,1)
a(k_shape.is_k_reducible_by_rectangle(s, k, (w,h)), True)

s = Partition([2, 1])
k = 1
(w,h) = (1,1)
a(k_shape.is_k_reducible_by_rectangle(s, k, (w,h)), True)

s = Partition([1, 1])
k = 2
(w,h) = (1,1)
a(k_shape.is_k_reducible_by_rectangle(s, k, (w,h)), False)

s = Partition([1, 1])
k = 2
(w,h) = (1,2)
a(k_shape.is_k_reducible_by_rectangle(s, k, (w,h)), True)

s = Partition([1, 1])
k = 2
(w,h) = (2,1)
a(k_shape.is_k_reducible_by_rectangle(s, k, (w,h)), False)

s = Partition([2, 1, 1])
k = 2
(w,h) = (2,1)
a(k_shape.is_k_reducible_by_rectangle(s, k, (w,h)), False)

s = Partition([2, 1, 1])
k = 2
(w,h) = (1,2)
a(k_shape.is_k_reducible_by_rectangle(s, k, (w,h)), True)

s = Partition([2, 1, 1])
k = 2
(w,h) = (1,1)
a(k_shape.is_k_reducible_by_rectangle(s, k, (w,h)), True)

s = Partition([2, 1, 1])
k = 3
(w,h) = (3,1)
a(k_shape.is_k_reducible_by_rectangle(s, k, (w,h)), False)

s = Partition([2, 1, 1])
k = 3
(w,h) = (2,2)
a(k_shape.is_k_reducible_by_rectangle(s, k, (w,h)), False)

s = Partition([2, 1, 1])
k = 3
(w,h) = (1,3)
a(k_shape.is_k_reducible_by_rectangle(s, k, (w,h)), False)

s = Partition([2, 1, 1])
k = 3
(w,h) = (2,1)
a(k_shape.is_k_reducible_by_rectangle(s, k, (w,h)), False)

s = Partition([2, 1, 1])
k = 3
(w,h) = (1,2)
a(k_shape.is_k_reducible_by_rectangle(s, k, (w,h)), True)

s = Partition([3, 2, 1])
k = 3
(w,h) = (3,1)
a(k_shape.is_k_reducible_by_rectangle(s, k, (w,h)), False)

s = Partition([3, 2, 1])
k = 3
(w,h) = (2,2)
a(k_shape.is_k_reducible_by_rectangle(s, k, (w,h)), True)

s = Partition([3, 2, 1])
k = 3
(w,h) = (1,3)
a(k_shape.is_k_reducible_by_rectangle(s, k, (w,h)), False)

s = Partition([3, 2, 1])
k = 3
(w,h) = (2,1)
a(k_shape.is_k_reducible_by_rectangle(s, k, (w,h)), False)

s = Partition([3, 2, 1])
k = 3
(w,h) = (1,2)
a(k_shape.is_k_reducible_by_rectangle(s, k, (w,h)), False)


# test_k_shape.is_reducible
s = Partition([1])
k = 1
a(k_shape.is_reducible(s, k), True)

s = Partition([2, 1])
k = 1
a(k_shape.is_reducible(s, k), True)

s = Partition([1, 1])
k = 2
a(k_shape.is_reducible(s, k), True)

s = Partition([2, 1, 1])
k = 2
a(k_shape.is_reducible(s, k), True)

s = Partition([2, 1, 1])
k = 3
a(k_shape.is_reducible(s, k), True)

s = Partition([3, 2, 1])
k = 3
a(k_shape.is_reducible(s, k), True)

s = Partition([5, 3, 2, 1, 1])
k = 4
a(k_shape.is_reducible(s, k), False)

s = Partition([5, 4, 2, 2, 1])
k = 4
a(k_shape.is_reducible(s, k), False)


# test_k_to_irreducible_k_shapes
ptns = k_to_irreducible_k_shapes(1)
a(ptns, [[]])

ptns = k_to_irreducible_k_shapes(2)
a(ptns, [[]])

ptns = k_to_irreducible_k_shapes(3)
a(ptns, [[], [1], [2, 1]])


# test is k core, now using builtin is core
a(Partition([2, 1]).is_core(1), False)
a(Partition([2, 1]).is_core(2), True)
a(Partition([2, 1]).is_core(3), False)
a(Partition([2, 1]).is_core(4), True)
a(Partition([4, 2, 2]).is_core(1), False)
a(Partition([4, 2, 2]).is_core(2), False)
a(Partition([4, 2, 2]).is_core(3), False)
a(Partition([4, 2, 2]).is_core(4), True)
a(Partition([4, 2, 2]).is_core(5), False)
a(Partition([4, 2, 2]).is_core(6), False)
a(Partition([4, 2, 2]).is_core(7), True)


# test partition to k core (or k bounded partition to k+1 core)
# main tests
a(to_k_core([1], 3), [1])
a(to_k_core([2], 3), [2])
a(to_k_core([1, 1], 3), [1, 1])
a(to_k_core([2, 1], 3), [3, 1])
a(to_k_core([1, 1, 1], 3), [2, 1, 1])
a(to_k_core([1, 1, 1, 1, 1, 1], 2), [6, 5, 4, 3, 2, 1])
a(to_k_core([2, 1, 1, 1], 3), [4, 2, 1, 1])
a(to_k_core([1, 1, 1, 1, 1], 3), [3, 2, 2, 1, 1])
# extra
a(to_k_core([6, 5, 5, 2], 4), [11, 8, 5, 2])
a(to_k_core([2, 2, 1], 3), [5, 3, 1])


# test_is_linked
# empty skew
sp = SkewPartition([[], []])
a(is_linked(sp), True)

# valid row shape but invalid col shape
sp = SkewPartition([[3, 2], [1, 0]])
a(is_linked(sp), False)

###########################################################################
# test_add_row
# sp = SkewPartition([[4, 3, 1], [1, 1]])
# a(add_row(sp, 1, 2), [[6, 5, 3, 1], [3, 3, 2]])


# test_thing_to_added_row_things
# sp = SkewPartition([[1], []])
# a(thing_to_added_row_things(sp, 0), [[[1],[]]])

# sp = SkewPartition([[5, 2, 1], [2]])
# a(thing_to_added_row_things(sp, 0), [[[5, 2, 1], [2]]])

# sp = SkewPartition([[4, 3, 1], [1, 1]])
# a(thing_to_added_row_things(sp, 0), [])

# sp = SkewPartition([[1], []])
# a(thing_to_added_row_things(sp, 1), [[[1, 1],[]], [[2, 1],[1]]])

# sp = SkewPartition([[2], []])
# a(thing_to_added_row_things(sp, 2), [[[2, 2],[]], [[3, 2],[1]], [[4, 2],[2]]])


# test_ptn_to_linked_things
# p = []
# a(ptn_to_linked_things(p), [[[],[]]])

# p = [7]
# a(ptn_to_linked_things(p), [[[7],[]]])

# p = [1, 1]
# a(ptn_to_linked_things(p), [[[1, 1],[]], [[2, 1],[1]]])

# p = [3, 2, 1]
# a(ptn_to_linked_things(p), [[[3, 2, 1],[]], [[4, 3, 1],[1, 1]], [[4, 2, 1],[1]], [[5, 2, 1], [2]], [[6, 3, 1],[3, 1]]])


# test_row_shape_to_linked_skew_partitions
p = []
a(row_shape_to_linked_skew_partitions(p), [[[],[]]])

p = [7]
a(row_shape_to_linked_skew_partitions(p), [[[7],[]]])

p = [1, 1]
a(row_shape_to_linked_skew_partitions(p), [[[1, 1],[]], [[2, 1],[1]]])

p = [3, 2, 1]
a(row_shape_to_linked_skew_partitions(p), [[[3, 2, 1],[]], [[4, 2, 1],[1]], [[5, 2, 1], [2]], [[6, 3, 1],[3, 1]]])

p = [2, 2, 2]
a(row_shape_to_linked_skew_partitions(p), [[[2, 2, 2],[]], [[4, 2, 2],[2]], [[6, 4, 2],[4, 2]]])

p = [3, 1, 1]
a(row_shape_to_linked_skew_partitions(p), [[[3, 1, 1],[]], [[4, 1, 1],[1]], [[5, 2, 1],[2, 1]]])


# test_complement
ri = RootIdeal([], n=1)
a(ri.complement(), [])

ri = RootIdeal([], n=2)
a(ri.complement(), [(0,1)])
a(ri.complement().complement(), [])

ri = RootIdeal([(0, 1)])
a(ri.complement().n, 2)

ri = RootIdeal([(0,1)], n=2)
a(ri.complement(), [])

ri = RootIdeal([(0,2), (0,3), (0,4), (1,3), (1,4), (2,3), (2,4)], n=5)
a(ri.complement(), [(0,1), (1,2), (3,4)])

ri = RootIdeal([(0,1), (1,2), (3,4)], n=5)
a(ri.complement(), [(0,2), (0,3), (0,4), (1,3), (1,4), (2,3), (2,4)])


# test_partition_to_k_schur_root_ideal
p = [2, 1]
n = 4
k = 2
a(RIS.init_k_schur_from_pseudo_partition(p, k, n), [(0,1), (0,2), (0,3), (1,3)])

p = [2, 1]
n = 4
k = 3
a(RIS.init_k_schur_from_pseudo_partition(p, k, n), [(0,2), (0,3)])

p = [2, 1]
n = 4
k = 4
a(RIS.init_k_schur_from_pseudo_partition(p, k, n), [(0,3)])

p = [2, 1]
n = 4
k = 5
a(RIS.init_k_schur_from_pseudo_partition(p, k, n), [])


# test compositional hall littlewood polynomial
Sym = SymmetricFunctions(QQ['t'])
hl = Sym.hall_littlewood().Qp()
for lis in ([3, 3, 2], [0], [1], [2], [1, 1], [2, 1], [2, 2, 1], [2, 1, 1], [6, 4, 2], [2, 2, 1, 1], [5, 3, 1], [5, 5, 3, 3, 1, 1], [4, 3, 2], [4, 4, 2], [4, 4, 1], [3, 1, 1], [4, 2, 2], [5, 5, 2, 2], [3, 3, 1, 1], [0], [1, 0], [2, 0], [3, 0]):
	p = Partition(lis)
	a(compositional_hall_littlewood_Qp(p), hl(p))

t = Sym.base_ring().gen()
HLQp = Sym.hall_littlewood().Qp()
# m + 1 = n
a(compositional_hall_littlewood_Qp([0, 1]), t * HLQp[1])
a(compositional_hall_littlewood_Qp([1, 2]), t * HLQp[2, 1])
a(compositional_hall_littlewood_Qp([2, 3]), t * HLQp[3, 2])
# m + 2 = n
a(compositional_hall_littlewood_Qp([0, 2]), t * HLQp[1, 1] + t * HLQp[2] - HLQp[1, 1])
a(compositional_hall_littlewood_Qp([0, 3]), t**2 * HLQp[2, 1] + t * HLQp[3] - HLQp[2, 1])
a(compositional_hall_littlewood_Qp([1, 3]), t * HLQp[2, 2] + t * HLQp[3, 1] - HLQp[2, 2])
# m + 4 = n
a(compositional_hall_littlewood_Qp([0, 4]), t**2 * HLQp[2, 2] + t**2 * HLQp[3, 1] - t * HLQp[2, 2] + t * HLQp[4] - HLQp[3, 1])
# length 3 composition
a(compositional_hall_littlewood_Qp([1, 3, 2]), t * HLQp[2, 2, 2] + t**2 * HLQp[3, 2, 1] - HLQp[2, 2, 2])


# test_straighten
Sym = SymmetricFunctions(QQ)
s = Sym.s()
a(straighten(s, [2, 1, 3]), -s[2, 2, 2])
h = Sym.h()
a(straighten(h, [5, 1, 7]), h[7, 5, 1])
e = Sym.e()
a(straighten(e, [5, 1, 7]), e[7, 5, 1])
p = Sym.p()
a(straighten(p, [5, 1, 7]), p[7, 5, 1])
w = Sym.w()
a(straighten(w, [5, 1, 7]), w[7, 5, 1])
# please verify that HOP(hl[3,0,1]) = t * hl[3,1]
base_ring = QQ['t']
Sym = SymmetricFunctions(base_ring)
t = base_ring.gen()
hl = Sym.hall_littlewood().Qp()
a(straighten(hl, [3, 0, 1]), t*hl[3, 1])


# seq space
S = ShiftingSequenceSpace()
a((1, -1) in S, True)
a((1, -1, 0, 9) in S, True)
a((0, 1, -1, 0, 9) in S, True)
a((1, -1, 0, 9, -9) in S, True)
a((1, -1, 0, -9, 9) in S, True)
a((1, 2, 1, -1, -2, 0, -1) in S, True)
a([1, -1] in S, False)
a((0.5,) in S, False)


# raising sequence space
S = RaisingSequenceSpace()
a((1, -1) in S, True)
a((1, -1, 0, 9) in S, False)
a((1, -1, 0, 9, -9) in S, True)
a((1, -1, 0, -9, 9) in S, False)
a((1, 2, 1, -1, -2, 0, -1) in S, True)
a([1, -1] in S, False)
a((0.5,) in S, False)


# test ShiftingOperatorAlgebra
R = ShiftingOperatorAlgebra()
TestSuite(R).run(verbose=False)
a(R(), R())
a(R[(1, 2, 1)] * R[(0, 1, 0, 1)], R[(1, 3, 1, 1)])
# retrieve indices
a(R[(1, -1)].indices(), [(1, -1)])
a((R[(1, -1)] + 2 * R[(1, 0, -1)]).indices(), [(1, -1), (1, 0, -1)])
# retrieve index
a(R[(1,-1)].index(), (1,-1))
# act on lists
a(R[(1, 0, -3)]([2, 2]), [([3, 2, -3], 1)])
a(R[(2, 1, -1, -2)]([2, 2]), [([4, 3, -1, -2], 1)])
# act on symmetric functions
base_ring = QQ['t']
Sym = SymmetricFunctions(base_ring)
t = base_ring.gen()
# act on s
s = Sym.s()
a(R[(1, -1)](s[2, 1]), s[3])
a(R[(1, -1)](t * s[2, 1]), t * s[3])
# act on s with straightening
a(R[(1, 0, 2)](s[1, 1, 1]), -s[2, 2, 2])
# act on h
h = Sym.h()
a(R[(1, -1)](h[2, 1]), h[3])
# act on HL Q'
hl = Sym.hall_littlewood().Qp()
a(R[(1, -1)](hl[2, 1]), hl[3])
# R() - t*R(0, 1, -1) - t*R(1, -1) + (t^2-t)*R(1, 0, -1) + t^2*R(1, 1, -2) + t^2*R(2, -1, -1) - t^3*R(2, 0, -2)
a(R.one()(hl[2, 1, 1]), hl[2, 1, 1])
a((-t*R[(0, 1, -1)])(hl[2, 1, 1]), -t * hl[2, 2])
a((R[(1, -1)])(hl[2, 1, 1]), t * hl[3, 1])
a((t*R[(1, -1)])(hl[2, 1, 1]), t**2 * hl[3, 1])
a((-t*R[(1, -1)])(hl[2, 1, 1]), -t**2 * hl[3, 1])
a(((t**2-t)*R[(1, 0, -1)])(hl[2, 1, 1]), (t**2-t)*hl[3, 1])
a((t**2*R[(1, 1, -2)])(hl[2, 1, 1]), 0)
a((t**2*R[(2, -1, -1)])(hl[2, 1, 1]), t**2*hl[4])
a((t**3*R[(2, 0, -2)])(hl[2, 1, 1]), 0)
a((-t**3*R[(2, 0, -2)])(hl[2, 1, 1]), 0)
# act on tuple of lists
a(R[(1, 0, -3)](([2, 2], [1], [-1])), ([([3, 2, -3], 1)], [([2, 0, -3], 1)], [([0, 0, -3], 1)]))
# act on tuple of partitions
s = Sym.s()
a(R[(3, 2, 2)]((s[4, 4], s[1], s[2, 1])), (s[7, 6, 2], s[4, 2, 2], s[5, 3, 2]))
# act on things added together
s = Sym.s()
a(R[(1, -1)](s[2, 1] + s[3, 1]), s[3] + s[4])


# test RaisingOperatorAlgebra
R = RaisingOperatorAlgebra()
TestSuite(R).run(verbose=False)
a(R[(1, -1)] * R[(0, 1, 0, -1)], R[(1, 0, 0, -1)])
# create 'R_ij' element
a(R.ij(1, 3), R[(0, 1, 0, -1)])
# retrieve indices
a(R[(1, -1)].indices(), [(1, -1)])
a((R[(1, -1)] + 2 * R[(1, 0, -1)]).indices(), [(1, -1), (1, 0, -1)])
# retrieve index
a(R[(1, -1)].index(), (1, -1))
# R(R(h)) does NOT necessarily equal (R*R)(h)
Sym = SymmetricFunctions(QQ)
h = Sym.h()
a((R.ij(0, 1) * R.ij(0, 1))(h[1, 1, 1, 1]), 0)
a((R.ij(1, 2) * R.ij(1, 3) * R.ij(0, 1) * R.ij(0, 1))(h[1, 1, 1, 1]), h[3, 1])


# test infinite dimensional free algebra ('x' variables indexed by NN)
R = InfiniteDimensionalFreeAlgebra()
TestSuite(R).run(verbose=False)
a(R(), R())
a(R(), R(0))
a(R.one(), R(1))
a(2 * R.one(), R(2))
a(R[1] + R[1], R(2 * R[1]))
a(5 * R[3] + R[4], 5 * R[3] + R[4])
a(R[1] * R[2], R[2] * R[1])
a(R[1] + R[2], R[2] + R[1])
# retrieve indices
# not implemented
# retrieve index
# not implemented


# test E
# e = InfiniteDimensionalFreeRing(prefix='a', index_set=IntegerRing())
# print(e)


# test double ring
DR = DoubleRing
TestSuite(DR).run(verbose=False)
a(5 * DR[3] + DR[-4], DR[-4] + DR[3] * 5)


# test hall littlewood vertex operator
Sym = SymmetricFunctions(QQ['t'])
hl = Sym.hall_littlewood().Qp()
H = HallLittlewoodVertexOperator
(g1, g2, g3) = (4, 1, 3)
gamma = [g1, g2, g3]
one = hl.one()
a(H(gamma)(one), H(g1)(H(g2)(H(g3)(one))))


# test qt raising roots operator
base_ring = QQ['t', 'q']
(t, q) = base_ring.gens()
s = SymmetricFunctions(base_ring).s()
op = qt_raising_roots_operator(3)
a(op(s[4, 2]), s[4, 2] + (-t-q)*s[5, 1] + t*q*s[6])


# test catalan function init methods
CFS = CatalanFunctions()

cf = CFS.init_from_indexed_root_ideal([(0,2), (1,2)], [6, 6, 5])
a(cf.roots, [(0,2), (1,2)])
a(cf.index, [6, 6, 5])

K = FractionField(QQ['t'])
elm = CFS.init_from_indexed_root_ideal([], [3, 3, 2, 1], base_ring=K)


# test catalan function
base_ring = QQ['t']
t = base_ring.gen()
sym = SymmetricFunctions(base_ring)
hl = sym.hall_littlewood().Qp()
spec_hl = sym.hall_littlewood(t=1).Qp()
s = sym.s()
# empty product
ri = RIS.init_from_partition([2, 1], n=3)
g = [3, 1, 1]
cat_func = CatalanFunction(ri, g)
a(cat_func.eval(), hl[3, 1, 1])
a(cat_func.eval(t=1), spec_hl[3, 1, 1])

for gamma in [[1], [1, 1], [2], [2, 1]]:
        a(cf.eval(t=1), spec_hl(gamma))

gamma = [1]
cf = CatalanFunction([], gamma)
a(cf.eval(), hl(gamma))

gamma = [1, 1]
cf = CatalanFunction([(0,1)], gamma)
a(cf.eval(), hl(gamma))

gamma = [2]
cf = CatalanFunction([], gamma)
a(cf.eval(), hl(gamma))

gamma = [2, 1]
cf = CatalanFunction([(0,1)], gamma)
a(cf.eval(), hl(gamma))

# if Psi is empty, then H(Psi, gamma) = s_gamma
gammas = ([1], [1, 1], [2], [2, 1], [3, 1], [4, 1], [2, 2], [3, 3], [5, 5], [5, 4], [6, 4], [10, 2], [3, 3, 1], [3, 1, 1], [5, 4, 3], [4, 2, 1], [5, 3, 1], [3, 3, 2], [3, 2, 2], [4, 4, 2], [4, 2, 2], [2, 2, 1], [2, 1, 1], [2, 2, 1, 1])
for gamma in gammas:
	cf = CatalanFunction([], gamma)
	a(s(cf.eval()), s(gamma))

ri = RIS.init_from_partition([1, 1], n=3)
g = [3, 1, 1]
cat_func = CatalanFunction(ri, g)
a(cat_func.eval(), hl[3, 1, 1] - t**2*hl[4, 1])

# issue #9
K = FractionField(QQ['t'])
cf = CFS.init_from_indexed_root_ideal([], [3, 3, 2, 1], base_ring=K)
cf.eval()


# test catalan function expand
cf = CatalanFunction([], [4, 1])
a(cf.expand(1), 0)
R = PolynomialRing(QQ, 2, 'x')
x = R.gens()
a(cf.expand(2), x[0]**4*x[1] + x[0]**3*x[1]**2 + x[0]**2*x[1]**3 + x[0]*x[1]**4)


# test staircase shape
a(staircase_shape(3), [2, 1, 0])


# test staircase root ideal
a(RIS.init_staircase(3), [(0,1), (0,2), (1,2)])
a(RIS.init_staircase(4), [(0,1), (0,2), (0,3), (1,2), (1,3), (2,3)])


# test parabolic root ideal
a(RIS.init_parabolic_from_composition([1]), [])
a(RIS.init_parabolic_from_composition([1]).n, 1)
a(RIS.init_parabolic_from_composition([2]), [])
a(RIS.init_parabolic_from_composition([2]).n, 2)
a(RIS.init_parabolic_from_composition([1, 1]), [(0, 1)])
a(RIS.init_parabolic_from_composition([1, 2]), [(0, 1), (0, 2)])
a(RIS.init_parabolic_from_composition([2, 1]), [(0, 2), (1, 2)])
a(RIS.init_parabolic_from_composition([2, 2]), [(0, 2), (0, 3), (1, 2), (1, 3)])
a(RIS.init_parabolic_from_composition([1, 3, 2]), [(0,1), (0,2), (0,3), (0,4), (0,5), (1,4), (1,5), (2,4), (2,5), (3,4), (3,5)])
# and different input types
a(RIS.init_parabolic_from_composition([1, 1]), [(0, 1)])
a(RIS.init_parabolic_from_composition(Partition([1, 1])), [(0, 1)])
a(RIS.init_parabolic_from_composition(Composition([1, 1])), [(0, 1)])


# test bottom for root ideal
ri = RIS.init_from_partition([3, 2, 1], 5)
a(ri.bottom(0), 4)
a(ri.bottom(1), 3)
a(ri.bottom(2), 4)
a(ri.bottom(3), 3)
a(ri.bottom(4), 4)


# test down path
ri = RIS.init_from_partition([3, 2, 1], 5)
a(ri.down_path(0), [0, 2, 4])


# test dual k theoretic h
a(dual_k_theoretic_homogeneous(0, 0), 1)
a(dual_k_theoretic_homogeneous(0, 13), 1)

Sym = SymmetricFunctions(QQ)
h = Sym.h()
a(dual_k_theoretic_homogeneous(1, 1), h[1] + 1)

Sym = SymmetricFunctions(QQ['t'])
h = Sym.h()
a(dual_k_theoretic_homogeneous(1, 1, base_ring=QQ['t']), h[1] + 1)

Sym = SymmetricFunctions(QQ['t'])
h = Sym.h()
a(dual_k_theoretic_homogeneous(1, 2, base_ring=QQ['t']), h[1] + 2)

Sym = SymmetricFunctions(QQ['t'])
h = Sym.h()
a(dual_k_theoretic_homogeneous(2, 1, base_ring=QQ['t']), h[2] + h[1] + 1)

Sym = SymmetricFunctions(QQ['t'])
h = Sym.h()
a(dual_k_theoretic_homogeneous(2, 2, base_ring=QQ['t']), h[2] + 2*h[1] + 3)

Sym = SymmetricFunctions(QQ['t'])
h = Sym.h()
a(dual_k_theoretic_homogeneous(2, 3, base_ring=QQ['t']), h[2] + 3*h[1] + 6)

# h_[2,1](x, [1, 1])
Sym = SymmetricFunctions(QQ)
h = Sym.h()
a(dual_k_theoretic_homogeneous([2, 1], [1, 1]), h[1]**2 + h[1]*h[2] + 2*h[1] + h[2] + 1)

# h_[1, 2](x, [2, 3])
Sym = SymmetricFunctions(QQ)
h = Sym.h()
a(dual_k_theoretic_homogeneous([1, 2], [2, 3]), 3*h[1]**2 + h[1]*h[2] + 2*h[2] + 12*h[1] + 12)


# test dual grothendieck function
Sym = SymmetricFunctions(QQ)
h = Sym.h()
func = dual_grothendieck_function([2, 1])
a(func, h[1]*h[2] + h[2] - h[3])

s = Sym.s()
sfunc = s(func)

Sym = SymmetricFunctions(QQ['t'])
h = Sym.h()
func = dual_grothendieck_function([2, 1], base_ring=QQ['t'])
a(func, h[1]*h[2] + h[2] - h[3])

s = Sym.s()
sfunc = s(func)


# test pieri operator
u = PieriOperatorAlgebra()
a(u.i(0)([3, 2, 1]), [([2, 2, 1], 1)])
a(u.i(1)([3, 2, 1]), [([3, 1, 1], 1)])
# act on catalan function
cf = CatalanFunction([(0,2), (1,2)], [6, 6, 6])
a(u.i(2)(cf), CatalanFunction([(0,2), (1,2)], [6, 6, 5]))
# act on k schur function
base_ring = QQ['t']
Sym = SymmetricFunctions(base_ring)
t = base_ring.gen()
ks = Sym.kBoundedSubspace(4, t).kschur()
# TODO: verify by hand that below is really correct, or maybe a simpler example.
# a(u.i(2)(ks[2, 2, 1]), ks[2, 2, 1] + t**2*ks[3, 2] + t**3*ks[4, 1])


# test double h
# h = double_homogeneous_building_block(1, 1)
# print(h)


# test k coverees 1
a(k_coverees1([3, 3, 2, 2, 1, 1], 2), set([Partition([3, 2, 2, 1, 1])]))
a(k_coverees1([2, 2, 1, 1], 2), set([Partition([2, 1, 1])]))
a(k_coverees1([2, 1, 1], 2), set([Partition([2]), Partition([1, 1])]))

a(k_coverees1([6, 4, 2, 2, 1], 5), set([Partition([5, 4, 2, 2, 1]), Partition([6, 2, 2, 2, 1]), Partition([6, 3, 2, 2]), Partition([6, 4, 2, 1, 1])]))


# test k coverees
a(k_coverees([3, 3, 2, 2, 1, 1], 2), set([Partition([3, 2, 2, 1, 1])]))
a(k_coverees([2, 2, 1, 1], 2), set([Partition([2, 1, 1])]))
a(k_coverees([2, 1, 1], 2), set([Partition([2]), Partition([1, 1])]))

a(k_coverees([6, 4, 2, 2, 1], 5), set([Partition([5, 4, 2, 2, 1]), Partition([6, 2, 2, 2, 1]), Partition([6, 3, 2, 2]), Partition([6, 4, 2, 1, 1])]))


# test go to ribbon head
a(__go_to_ribbon_head([(0,1)], (0,1)), (0,1))
a(__go_to_ribbon_head([(1,0), (2,0)], (2,0)), (1,0))
a(__go_to_ribbon_head([(1,0), (2,0), (3,0)], (3,0)), (1,0))
a(__go_to_ribbon_head([(1,1), (1,2)], (1,1)), (1,2))
a(__go_to_ribbon_head([(1,1), (1,2), (1,3)], (1,1)), (1,3))
a(__go_to_ribbon_head([(1,1), (1,2), (2,1)], (2,1)), (1,2))
a(__go_to_ribbon_head([(0,3), (1,1), (1,2), (2,0), (2,1), (3,0), (3,1)], (2,1)), (1,2))


# test row marking to marking
a(row_marking_to_marking([3, 2, 2], [2, 1], 0), (0,2))
a(row_marking_to_marking([3, 2, 2], [2, 1], 1), (1,1))

a(row_marking_to_marking([3, 2, 2], [1, 1], 0), (0,2))

a(row_marking_to_marking([3, 3, 2, 2], [2, 2, 1, 1], 0), (0,2))
a(row_marking_to_marking([3, 3, 2, 2], [2, 2, 1, 1], 2), (2,1))


# test row markings to markings
a(row_markings_to_markings(([], [1]), [0]), [(0,0)])
a(row_markings_to_markings(([], [1], [1, 1]), [0, 1]), [(0,0), (1,0)])
a(row_markings_to_markings(([], [1], [1, 1], [2, 1, 1]), [0, 1, 2]), [(0,0), (1,0), (2,0)])
a(row_markings_to_markings(([], [1], [1, 1], [2, 1, 1], [3, 1, 1]), [0, 1, 2, 0]), [(0,0), (1,0), (2,0), (0,2)])
a(row_markings_to_markings(([], [1], [1, 1], [2, 1, 1], [3, 1, 1], [5, 3, 1]), [0, 1, 2, 0, 1]), [(0,0), (1,0), (2,0), (0,2), (1,2)])


# test is markable
a(is_row_markable([3, 2, 2], [2, 1], 0), True)
a(is_row_markable([3, 2, 2], [2, 1], 1), True)
a(is_row_markable([3, 2, 2], [2, 1], 2), False)
a(is_row_markable([3, 2, 2], [2, 1], 3), False)

a(is_row_markable([3, 2, 2], [1, 1], 0), True)
a(is_row_markable([3, 2, 2], [1, 1], 1), False)
a(is_row_markable([3, 2, 2], [1, 1], 2), False)
a(is_row_markable([3, 2, 2], [1, 1], 3), False)

a(is_row_markable([3, 3, 2, 2], [2, 2, 1, 1], 0), True)
a(is_row_markable([3, 3, 2, 2], [2, 2, 1, 1], 1), False)
a(is_row_markable([3, 3, 2, 2], [2, 2, 1, 1], 2), True)
a(is_row_markable([3, 3, 2, 2], [2, 2, 1, 1], 3), False)


# test k marked coverees
# a(k_coverees([6, 4, 2, 2, 1], 5), set([Partition([5, 4, 2, 2, 1]), Partition([6, 2, 2, 2, 1]), Partition([6, 3, 2, 2]), Partition([6, 4, 2, 1, 1])]))
a(k_marked_coverees([6, 4, 2, 2, 1], 5, 0), set([Partition([5, 4, 2, 2, 1])]))
a(k_marked_coverees([6, 4, 2, 2, 1], 5, 1), set([Partition([6, 2, 2, 2, 1]), Partition([6, 3, 2, 2])]))
a(k_marked_coverees([6, 4, 2, 2, 1], 5, 2), set())
a(k_marked_coverees([6, 4, 2, 2, 1], 5, 3), set([Partition([6, 4, 2, 1, 1])]))
a(k_marked_coverees([6, 4, 2, 2, 1], 5, 4), set([Partition([6, 3, 2, 2])]))


# test end core to marked core sequences
a(end_core_to_marked_core_sequences([5, 3, 1], 2, [0, 1, 2, 0, 1]),
	set([partition_tuple([], [1], [1, 1], [2, 1, 1], [3, 1, 1], [5, 3, 1])]))
a(end_core_to_marked_core_sequences([5, 3, 1], 2, [1, 2, 0, 1]),
	set([partition_tuple([1], [1, 1], [2, 1, 1], [3, 1, 1], [5, 3, 1])]))
a(end_core_to_marked_core_sequences([5, 3, 1], 2, [1]),
	set([
		partition_tuple([3, 1, 1], [5, 3, 1]),
		partition_tuple([4, 2], [5, 3, 1]),
	]))
a(end_core_to_marked_core_sequences([5, 3, 1], 2, [0]),
	set([
		partition_tuple([3, 1, 1], [5, 3, 1]),
		partition_tuple([4, 2], [5, 3, 1]),
	]))
a(end_core_to_marked_core_sequences([5, 3, 1], 2, [2]),
	set([
		partition_tuple([4, 2], [5, 3, 1]),
	]))
a(end_core_to_marked_core_sequences([5, 3, 1], 2, [0, 1]),
	set([
		partition_tuple([2, 1, 1], [3, 1, 1], [5, 3, 1]),
		partition_tuple([3, 1], [4, 2], [5, 3, 1]),
	]))
a(end_core_to_marked_core_sequences([2, 1, 1], 2, [2]),
	set([
		partition_tuple([1, 1], [2, 1, 1]),
	]))
a(end_core_to_marked_core_sequences([5, 3, 1], 2, [0, 2, 0]),
	set([
		partition_tuple([1, 1], [3, 1], [3, 1, 1], [5, 3, 1]),
		partition_tuple([2], [3, 1], [3, 1, 1], [5, 3, 1]),
	]))
a(end_core_to_marked_core_sequences([5, 3, 1], 2, [0, 1, 1, 2]),
	set([
		partition_tuple([1], [2], [3, 1], [4, 2], [5, 3, 1]),
	]))
a(end_core_to_marked_core_sequences([5, 3, 1], 2, [1, 2, 2]), set())


# test std_strong_tab_from_core_sequence(core_sequence, k, marks):
cs = [[], [1]]
k = 1
marks = [(0,0)]
st = std_strong_tab_from_core_sequence(cs, k, marks)
a(st, StrongTableau([[-1]], k))

cs = [[], [1], [2]]
k = 4
marks = [(0,0), (0,1)]
st = std_strong_tab_from_core_sequence(cs, k, marks)
a(st, StrongTableau([[-1, -2]], k))

cs = [[], [1], [2], [2, 1], [3, 1]]
k = 4
marks = [(0,0), (0,1), (1,0), (0,2)]
st = std_strong_tab_from_core_sequence(cs, k, marks)
a(st, StrongTableau([[-1, -2, -4], [-3]], k))


# test end core to strong marked tableaux
a(end_core_to_strong_marked_tableaux([5, 3, 1], 2, [0, 1, 2, 0, 1]),
	set([StrongTableau([[-1, 3, -4, 5, 5], [-2, 5, -5], [-3]], 2)]))
a(end_core_to_strong_marked_tableaux([5, 3, 1], 2, [1]),
	set([
		StrongTableau([[None, None, None, 1, 1], [None, 1, -1], [None]], 2),
		StrongTableau([[None, None, None, None, 1], [None, None, -1], [1]], 2),
	]))


# test ungraded
# setup
base_ring = ZZ['t']
t = base_ring.gen()
sym = SymmetricFunctions(base_ring)
s = sym.s()
h = sym.h()
hl = sym.hall_littlewood().Qp()

# test
f = t * s[2, 1]
a(ungraded(f), s[2, 1])

f = 1 + t**2 * s[1, 1] - 4 * t * s[2, 2, 1, 1]
a(ungraded(f), 1 + s[1, 1] - 4 * s[2, 2, 1, 1])

f = t * hl[2, 1]
a(ungraded(f), hl[2, 1])

f = 1 + t**2 * hl[1, 1] - 4 * t * hl[2, 2, 1, 1]
a(ungraded(f), 1 + hl[1, 1] - 4 * hl[2, 2, 1, 1])

f = t * h[2, 1]
a(ungraded(f), h[2, 1])

f = 1 + t**2 * h[1, 1] - 4 * t * h[2, 2, 1, 1]
a(ungraded(f), 1 + h[1, 1] - 4 * h[2, 2, 1, 1])








# ALL DONE!
print('Testing completed successfully!', end='')
end_time = time.time()
elapsed_time = end_time - start_time
print(' Elapsed time = {}'.format(elapsed_time))
