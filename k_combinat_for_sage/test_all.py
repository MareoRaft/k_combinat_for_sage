#!/usr/bin/env sage
# A place to test my functions
from __future__ import print_function
import time

from sage.all import *
from testing import *
from all import *
start_time = time.time()
print('Sage loaded.  Testing...')



# test_SP_is_symmetric
sp = SkewPartition([[], []])
b = SP.is_symmetric(sp)
a(b, True)

sp = SkewPartition([[1], []])
b = SP.is_symmetric(sp)
a(b, True)

sp = SkewPartition([[4,3,3,1], [1]])
b = SP.is_symmetric(sp)
a(b, True)

sp = SkewPartition([[4,3,3,1], [1,1]])
b = SP.is_symmetric(sp)
a(b, False)

sp = SkewPartition([[5,3,3,1], [2,2]])
b = SP.is_symmetric(sp)
a(b, False)



# test_right
sp = SkewPartition([[1], []])
c = SP.right(sp, 0)
a(c, 0)

sp = SkewPartition([[2], []])
c = SP.right(sp, 0)
a(c, 1)

sp = SkewPartition([[1, 1], []])
c = SP.right(sp, 0)
a(c, 0)

sp = SkewPartition([[2], [1]])
c = SP.right(sp, 0)
a(c, 1)

sp = SkewPartition([[2], [2]])
c = SP.right(sp, 0)
a(c, None)

sp = SkewPartition([[5, 5, 4, 2, 2],  [4, 3, 2]])
c = SP.right(sp, 2)
a(c, 3)


# test_left
sp = SkewPartition([[1], []])
c = SP.left(sp, 0)
a(c, 0)

sp = SkewPartition([[2], []])
c = SP.left(sp, 0)
a(c, 0)

sp = SkewPartition([[2], [1]])
c = SP.left(sp, 0)
a(c, 1)

sp = SkewPartition([[2], [2]])
c = SP.left(sp, 0)
a(c, None)

sp = SkewPartition([[5, 5, 4, 2, 2],  [4, 3, 2]])
c = SP.left(sp, 2)
a(c, 2)


# test_top
sp = SkewPartition([[1], []])
c = SP.top(sp, 0)
a(c, 0)

sp = SkewPartition([[2], []])
c = SP.top(sp, 0)
a(c, 0)

sp = SkewPartition([[1, 1], []])
c = SP.top(sp, 0)
a(c, 1)

sp = SkewPartition([[1, 1], [1]])
c = SP.top(sp, 0)
a(c, 1)

sp = SkewPartition([[2], [2]])
c = SP.top(sp, 0)
a(c, None)

sp = SkewPartition([[5, 5, 4, 2, 2],  [4, 3, 2]])
c = SP.top(sp, 2)
a(c, 2)

sp = SkewPartition([[5, 5, 4, 2, 2],  [4, 3, 2]])
c = SP.top(sp, 3)
a(c, 2)


# test_bottom
sp = SkewPartition([[1], []])
c = SP.bottom(sp, 0)
a(c, 0)

sp = SkewPartition([[2], []])
c = SP.bottom(sp, 0)
a(c, 0)

sp = SkewPartition([[1, 1], []])
c = SP.bottom(sp, 0)
a(c, 0)

sp = SkewPartition([[1, 1], [1]])
c = SP.bottom(sp, 0)
a(c, 1)

sp = SkewPartition([[2], [2]])
c = SP.bottom(sp, 0)
a(c, None)

sp = SkewPartition([[5, 5, 4, 2, 2],  [4, 3, 2]])
c = SP.bottom(sp, 2)
a(c, 2)

sp = SkewPartition([[5, 5, 4, 2, 2],  [4, 3, 2]])
c = SP.bottom(sp, 3)
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


# test_skew_partition_to_removable_roots
sp = SkewPartition([[], []])
a(skew_partition_to_removable_roots(sp, type='max'), [])

sp = SkewPartition([[3, 2, 1], [1]])
a(skew_partition_to_removable_roots(sp, type='max'), [(0,2)])

sp = SkewPartition([[6, 5, 3, 2, 2, 1], [2, 2]])
a(skew_partition_to_removable_roots(sp, type='max'), [(0,3), (1,4)])



# test_selected_rows_to_maximum_root_ideal
root_ideal = selected_rows_to_maximum_root_ideal(1, [0])
a(root_ideal, [])

root_ideal = selected_rows_to_maximum_root_ideal(2, [0, 1])
a(root_ideal, [])

root_ideal = selected_rows_to_maximum_root_ideal(2, [0])
a(root_ideal, [(0,1)])

root_ideal = selected_rows_to_maximum_root_ideal(3, [0, 2])
a(root_ideal, [(0,1), (0,2)])

root_ideal = selected_rows_to_maximum_root_ideal(4, [0, 2])
a(root_ideal, [(0,1), (0,2), (0,3), (1,3)])

root_ideal = selected_rows_to_maximum_root_ideal(5, [0, 1, 4])
a(root_ideal, [(0,2), (0,3), (0,4), (1,3), (1,4)])

root_ideal = selected_rows_to_maximum_root_ideal(5, [0, 1])
a(root_ideal, [(0,2), (0,3), (0,4), (1,3), (1,4), (2,4)])


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


# test_skew_partition_to_root_ideal
sp = SkewPartition([[6, 5, 3, 2, 2, 1], [2, 2]])
ri = skew_partition_to_root_ideal(sp, type='max')
a(ri, [(0,3), (0,4), (0,5), (1,4), (1,5)])

sp = SkewPartition([[4, 3, 2, 2, 1, 1], [3, 2, 1, 1]])
ri = skew_partition_to_root_ideal(sp, type='max')
a(ri, [(0,1), (0,2), (0,3), (0,4), (0,5), (1,2), (1,3), (1,4), (1,5), (2,4), (2,5), (3,5)])

sp = SkewPartition([[4, 3, 2, 2, 1, 1], [3, 2, 1, 1]])
ri = skew_partition_to_root_ideal(sp, type='min')
a(ri, [(0,1), (0,2), (0,3), (0,4), (0,5), (1,3), (1,4), (1,5), (2,4), (2,5), (3,5)])


# test_removable_roots_to_root_ideal
rr = [(0,1), (2,2)]
n = 4
a(removable_roots_to_root_ideal(rr, n), [(0,1), (0,2), (0,3), (1,2), (1,3), (2,2), (2,3)])


# test_down
# Rightmost example 2.4 of SKEW-LINKED CATALAN FUNCTIONS AND k-SCHUR POSITIVITY
root_ideal = [(0,1), (0,2), (0,3), (0,4), (0,5), (1,4), (1,5), (2,4), (2,5), (3,4), (3,5)]
a(down(root_ideal, 0), 1)
a(down(root_ideal, 1), 4)
a(down(root_ideal, 2), 4)
a(down(root_ideal, 3), 4)
a(down(root_ideal, 4), None)
a(down(root_ideal, 5), None)


# test_down_path
# Rightmost example 2.4 of SKEW-LINKED CATALAN FUNCTIONS AND k-SCHUR POSITIVITY
root_ideal = [(0,1), (0,2), (0,3), (0,4), (0,5), (1,4), (1,5), (2,4), (2,5), (3,4), (3,5)]
a(down_path(root_ideal, 0), [0, 1, 4])
a(down_path(root_ideal, 1), [1, 4])
a(down_path(root_ideal, 2), [2, 4])
a(down_path(root_ideal, 3), [3, 4])
a(down_path(root_ideal, 4), [4])
a(down_path(root_ideal, 5), [5])


# test_down_path_column_lengths_part
# Rightmost example 2.4 of SKEW-LINKED CATALAN FUNCTIONS AND k-SCHUR POSITIVITY
root_ideal = [(0,1), (0,2), (0,3), (0,4), (0,5), (1,4), (1,5), (2,4), (2,5), (3,4), (3,5)]
ptn = [7, 6, 5, 2, 2, 2]
a(down_path_column_lengths_part(root_ideal, ptn, 0), 15)
a(down_path_column_lengths_part(root_ideal, ptn, 1), 8)
a(down_path_column_lengths_part(root_ideal, ptn, 2), 7)
a(down_path_column_lengths_part(root_ideal, ptn, 3), 4)
a(down_path_column_lengths_part(root_ideal, ptn, 4), 2)
a(down_path_column_lengths_part(root_ideal, ptn, 5), 2)


# test_down_path_column_lengths
# Rightmost example 2.4 of SKEW-LINKED CATALAN FUNCTIONS AND k-SCHUR POSITIVITY
root_ideal = [(0,1), (0,2), (0,3), (0,4), (0,5), (1,4), (1,5), (2,4), (2,5), (3,4), (3,5)]
ptn = [7, 6, 5, 2, 2, 2]
a(down_path_column_lengths(root_ideal, ptn), [15, 7, 4, 2])


# test_root_ideal_to_partition
ri = []
a(root_ideal_to_partition(ri), [])

# Rightmost example 2.4 of SKEW-LINKED CATALAN FUNCTIONS AND k-SCHUR POSITIVITY
ri = [(0,1), (0,2), (0,3), (0,4), (0,5), (1,4), (1,5), (2,4), (2,5), (3,4), (3,5)]
a(root_ideal_to_partition(ri), [5, 2, 2, 2])


# test_partition_to_root_ideal
p = Partition([])
n = 5
a(partition_to_root_ideal(p, n), [])

# Rightmost example 2.4 of SKEW-LINKED CATALAN FUNCTIONS AND k-SCHUR POSITIVITY
p = Partition([5, 2, 2, 2])
n = 6
a(partition_to_root_ideal(p, n), [(0,1), (0,2), (0,3), (0,4), (0,5), (1,4), (1,5), (2,4), (2,5), (3,4), (3,5)])


# test_is_strict
ri = []
a(is_strict(ri), True)

# Rightmost example 2.4 of SKEW-LINKED CATALAN FUNCTIONS AND k-SCHUR POSITIVITY
ri = [(0,1), (0,2), (0,3), (0,4), (0,5), (1,4), (1,5), (2,4), (2,5), (3,4), (3,5)]
a(is_strict(ri), False)

ri = [(0,1), (0,2), (0,3), (0,4), (0,5)]
a(is_strict(ri), True)

ri = [(0,1), (0,2), (0,3), (0,4), (0,5), (1,2), (1,3), (1,4), (1,5), (2,3), (2,4), (2,5), (3,4), (3,5)]
a(is_strict(ri), True)

ri = [(0,1), (0,2), (0,3), (0,4), (0,5), (1,2), (1,3), (1,4), (1,5), (2,3), (2,4), (2,5), (3,4), (3,5), (4,5)]
a(is_strict(ri), True)


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


# test_P.next
min_ = []
max_ = []
p = []
a(P.next(p, min=min_, max=max_), func=is_False_or_None)

min_ = []
max_ = [1]
a(P.next([], min=min_, max=max_), [1])
a(P.next([1], min=min_, max=max_), func=is_False_or_None)

min_ = []
max_ = [1, 1]
a(P.next([], min=min_, max=max_), [1])
a(P.next([1], min=min_, max=max_), [1, 1])
a(P.next([1, 1], min=min_, max=max_), func=is_False_or_None)

min_ = []
max_ = [2]
a(P.next([], min=min_, max=max_), [1])
a(P.next([1], min=min_, max=max_), [2])
a(P.next([2], min=min_, max=max_), func=is_False_or_None)

min_ = []
max_ = [2, 1]
a(P.next([], min=min_, max=max_), [1])
a(P.next([1], min=min_, max=max_), [1, 1])
a(P.next([1, 1], min=min_, max=max_), [2])
a(P.next([2], min=min_, max=max_), [2, 1])
a(P.next([2, 1], min=min_, max=max_), func=is_False_or_None)

min_ = [1, 1]
max_ = [3, 2, 1]
a(P.next([1, 1], min=min_, max=max_), [1, 1, 1])
a(P.next([1, 1, 1], min=min_, max=max_), [2, 1])


# test_RootIdeal_next
ri = []
n = 4
a(RootIdeal_next(ri, n=n), [(0,3)])

max_ri = [(0,0), (0,1), (1,0), (1,1)]
min_ri = []
ri = [(0,1), (1,1)]
a(RootIdeal_next(ri, n=2, min=min_ri, max=max_ri), [(0,0), (0,1)])

max_ri = [(0,1), (0,2), (0,3), (1,2), (1,3), (2,3)]
min_ri = [(0,2), (0,3)]
ri = [(0,1), (0,2), (0,3)]
a(RootIdeal_next(ri, n=4, min=min_ri, max=max_ri), [(0,1), (0,2), (0,3), (1,3)])


# test_skew_partition_to_root_ideals
sp = SkewPartition([[4, 2, 1, 1], [2, 1]])
correct_ris = [
	RootIdeal([(0,1), (0,2), (0,3), (1,2), (1,3)]),
	RootIdeal([(0,1), (0,2), (0,3), (1,3)]),
]
a(set(skew_partition_to_root_ideals(sp)), set(correct_ris))

sp = SkewPartition([[4, 3, 2, 2, 1, 1], [3, 2, 1, 1]])
correct_ris = [
	RootIdeal([(0,1), (0,2), (0,3), (0,4), (0,5), (1,2), (1,3), (1,4), (1,5), (2,4), (2,5), (3,5)]),
	RootIdeal([(0,1), (0,2), (0,3), (0,4), (0,5), (1,3), (1,4), (1,5), (2,4), (2,5), (3,5)]),
]
a(set(skew_partition_to_root_ideals(sp)), set(correct_ris))


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


# test_kS.is_k_reducible_by_rectangle
s = Partition([1])
k = 1
(w,h) = (1,1)
a(kS.is_k_reducible_by_rectangle(s, k, (w,h)), True)

s = Partition([2, 1])
k = 1
(w,h) = (1,1)
a(kS.is_k_reducible_by_rectangle(s, k, (w,h)), True)

s = Partition([1, 1])
k = 2
(w,h) = (1,1)
a(kS.is_k_reducible_by_rectangle(s, k, (w,h)), False)

s = Partition([1, 1])
k = 2
(w,h) = (1,2)
a(kS.is_k_reducible_by_rectangle(s, k, (w,h)), True)

s = Partition([1, 1])
k = 2
(w,h) = (2,1)
a(kS.is_k_reducible_by_rectangle(s, k, (w,h)), False)

s = Partition([2, 1, 1])
k = 2
(w,h) = (2,1)
a(kS.is_k_reducible_by_rectangle(s, k, (w,h)), False)

s = Partition([2, 1, 1])
k = 2
(w,h) = (1,2)
a(kS.is_k_reducible_by_rectangle(s, k, (w,h)), True)

s = Partition([2, 1, 1])
k = 2
(w,h) = (1,1)
a(kS.is_k_reducible_by_rectangle(s, k, (w,h)), True)

s = Partition([2, 1, 1])
k = 3
(w,h) = (3,1)
a(kS.is_k_reducible_by_rectangle(s, k, (w,h)), False)

s = Partition([2, 1, 1])
k = 3
(w,h) = (2,2)
a(kS.is_k_reducible_by_rectangle(s, k, (w,h)), False)

s = Partition([2, 1, 1])
k = 3
(w,h) = (1,3)
a(kS.is_k_reducible_by_rectangle(s, k, (w,h)), False)

s = Partition([2, 1, 1])
k = 3
(w,h) = (2,1)
a(kS.is_k_reducible_by_rectangle(s, k, (w,h)), False)

s = Partition([2, 1, 1])
k = 3
(w,h) = (1,2)
a(kS.is_k_reducible_by_rectangle(s, k, (w,h)), True)

s = Partition([3, 2, 1])
k = 3
(w,h) = (3,1)
a(kS.is_k_reducible_by_rectangle(s, k, (w,h)), False)

s = Partition([3, 2, 1])
k = 3
(w,h) = (2,2)
a(kS.is_k_reducible_by_rectangle(s, k, (w,h)), True)

s = Partition([3, 2, 1])
k = 3
(w,h) = (1,3)
a(kS.is_k_reducible_by_rectangle(s, k, (w,h)), False)

s = Partition([3, 2, 1])
k = 3
(w,h) = (2,1)
a(kS.is_k_reducible_by_rectangle(s, k, (w,h)), False)

s = Partition([3, 2, 1])
k = 3
(w,h) = (1,2)
a(kS.is_k_reducible_by_rectangle(s, k, (w,h)), False)


# test_kS.is_reducible2
s = Partition([1])
k = 1
a(kS.is_reducible2(s, k), True)

s = Partition([2, 1])
k = 1
a(kS.is_reducible2(s, k), True)

s = Partition([1, 1])
k = 2
a(kS.is_reducible2(s, k), True)

s = Partition([2, 1, 1])
k = 2
a(kS.is_reducible2(s, k), True)

s = Partition([2, 1, 1])
k = 3
a(kS.is_reducible2(s, k), True)

s = Partition([3, 2, 1])
k = 3
a(kS.is_reducible2(s, k), True)

s = Partition([5, 3, 2, 1, 1])
k = 4
a(kS.is_reducible2(s, k), False)

s = Partition([5, 4, 2, 2, 1])
k = 4
a(kS.is_reducible2(s, k), False)


# test_k_to_irreducible_k_shapes
ptns = k_to_irreducible_k_shapes(1)
a(ptns, [[]])

ptns = k_to_irreducible_k_shapes(2)
a(ptns, [[]])

ptns = k_to_irreducible_k_shapes(3)
a(ptns, [[], [1], [2, 1]])


# test_is_k_core
a(is_k_core(Partition([2, 1]), 1), False)
a(is_k_core(Partition([2, 1]), 2), True)
a(is_k_core(Partition([2, 1]), 3), False)
a(is_k_core(Partition([2, 1]), 4), True)
a(is_k_core(Partition([4, 2, 2]), 1), False)
a(is_k_core(Partition([4, 2, 2]), 2), False)
a(is_k_core(Partition([4, 2, 2]), 3), False)
a(is_k_core(Partition([4, 2, 2]), 4), True)
a(is_k_core(Partition([4, 2, 2]), 5), False)
a(is_k_core(Partition([4, 2, 2]), 6), False)
a(is_k_core(Partition([4, 2, 2]), 7), True)


# test_is_linked
# empty skew
sp = SkewPartition([[], []])
a(is_linked(sp), True)

# valid row shape but invalid col shape
sp = SkewPartition([[3, 2], [1, 0]])
a(is_linked(sp), False)

###########################################################################
# test_add_row
sp = SkewPartition([[4, 3, 1], [1, 1]])
a(add_row(sp, 1, 2), [[6, 5, 3, 1], [3, 3, 2]])


# test_thing_to_added_row_things
sp = SkewPartition([[1], []])
a(thing_to_added_row_things(sp, 0), [[[1],[]]])

sp = SkewPartition([[5, 2, 1], [2]])
a(thing_to_added_row_things(sp, 0), [[[5, 2, 1], [2]]])

sp = SkewPartition([[4, 3, 1], [1, 1]])
a(thing_to_added_row_things(sp, 0), [])

sp = SkewPartition([[1], []])
a(thing_to_added_row_things(sp, 1), [[[1, 1],[]], [[2, 1],[1]]])

sp = SkewPartition([[2], []])
a(thing_to_added_row_things(sp, 2), [[[2, 2],[]], [[3, 2],[1]], [[4, 2],[2]]])


# test_ptn_to_linked_things
p = []
a(ptn_to_linked_things(p), [[[],[]]])

p = [7]
a(ptn_to_linked_things(p), [[[7],[]]])

p = [1, 1]
a(ptn_to_linked_things(p), [[[1, 1],[]], [[2, 1],[1]]])

p = [3, 2, 1]
a(ptn_to_linked_things(p), [[[3, 2, 1],[]], [[4, 3, 1],[1, 1]], [[4, 2, 1],[1]], [[5, 2, 1], [2]], [[6, 3, 1],[3, 1]]])


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
ri = []
n = 1
a(RI.complement(ri, n), [])

ri = []
n = 2
a(RI.complement(ri, n), [(0,1)])
a(RI.complement(RI.complement(ri, n), n), [])

ri = [(0,1)]
a(RI.complement(ri, n=2), [])

ri = [(0,2), (0,3), (0,4), (1,3), (1,4), (2,3), (2,4)]
a(RI.complement(ri, n=5), [(0,1), (1,2), (3,4)])

ri = [(0,1), (1,2), (3,4)]
a(RI.complement(ri, 5), [(0,2), (0,3), (0,4), (1,3), (1,4), (2,3), (2,4)])


# test_partition_to_k_schur_root_ideal
p = [2, 1]
n = 4
k = 2
a(partition_to_k_schur_root_ideal(p, k, n), [(0,1), (0,2), (0,3), (1,3)])

p = [2, 1]
n = 4
k = 3
a(partition_to_k_schur_root_ideal(p, k, n), [(0,2), (0,3)])

p = [2, 1]
n = 4
k = 4
a(partition_to_k_schur_root_ideal(p, k, n), [(0,3)])

p = [2, 1]
n = 4
k = 5
a(partition_to_k_schur_root_ideal(p, k, n), [])


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
Sym = SymmetricFunctions(QQ['t'])
# act on s
s = Sym.s()
a(R[(1, -1)](s[2, 1]), s[3])
# act on s with straightening
a(R[(1, 0, 2)](s[1, 1, 1]), -s[2, 2, 2])
# act on h
h = Sym.h()
a(R[(1, -1)](h[2, 1]), h[3])
# act on HL Q'
hl = Sym.hall_littlewood().Qp()
a(R[(1, -1)](hl[2, 1]), hl[3])
# act on tuple of lists
a(R[(1, 0, -3)](([2, 2], [1], [-1])), ([([3, 2, -3], 1)], [([2, 0, -3], 1)], [([0, 0, -3], 1)]))
# act on tuple of partitions
s = Sym.s()
a(R[(3, 2, 2)]((s[4, 4], s[1], s[2, 1])), (s[7, 6, 2], s[4, 2, 2], s[5, 3, 2]))


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


# test double ring
DR = DoubleRing
TestSuite(DR).run(verbose=False)
a(5 * DR[3] + DR[-4], DR[-4] + DR[3] * 5)


# test_straighten
Sym = SymmetricFunctions(QQ)
s = Sym.s()
a(straighten(s, [2, 1, 3]), -s[2, 2, 2])


# test hall littlewood vertex operator
Sym = SymmetricFunctions(QQ['t'])
hl = Sym.hall_littlewood().Qp()
H = HallLittlewoodVertexOperator
(g1, g2, g3) = (4, 1, 3)
gamma = [g1, g2, g3]
one = hl.one()
a(H(gamma)(one), H(g1)(H(g2)(H(g3)(one))))


# test compositional hall littlewood polynomial
Sym = SymmetricFunctions(QQ['t'])
hl = Sym.hall_littlewood().Qp()
a(compositional_hall_littlewood_Qp([3, 3, 2]), hl[3, 3, 2])


# test_indexed_root_ideal_to_catalan_function
Sym = SymmetricFunctions(QQ['t'])
hl = Sym.hall_littlewood().Qp()
# empty product
ri = partition_to_root_ideal([2, 1], n=3)
g = [3, 1, 1]
cat_func = indexed_root_ideal_to_catalan_function(ri, g)
a(cat_func, hl[3, 1, 1])
# other
# ri = partition_to_root_ideal([1, 1], n=3)
# g = [3, 1, 1]
# cat_func = indexed_root_ideal_to_catalan_function(ri, g)
# a(cat_func, hl[3, 1, 1])


# TODO: test all catalan function methods with some example.  Ask Jennifer morse.


# test staircase shape
a(staircase_shape(3), [2, 1, 0])


# test staircase root ideal
a(staircase_root_ideal(3), [(0,1), (0,2), (1,2)])
a(staircase_root_ideal(4), [(0,1), (0,2), (0,3), (1,2), (1,3), (2,3)])


# test dual k theoretic h
a(dual_k_theoretic_h(0, 0), 1)
a(dual_k_theoretic_h(0, 13), 1)

Sym = SymmetricFunctions(QQ)
h = Sym.h()
a(dual_k_theoretic_h(1, 1), h[1] + 1)

Sym = SymmetricFunctions(QQ['t'])
h = Sym.h()
a(dual_k_theoretic_h(1, 1, base_ring=QQ['t']), h[1] + 1)

Sym = SymmetricFunctions(QQ['t'])
h = Sym.h()
a(dual_k_theoretic_h(1, 2, base_ring=QQ['t']), h[1] + 2)

Sym = SymmetricFunctions(QQ['t'])
h = Sym.h()
a(dual_k_theoretic_h(2, 1, base_ring=QQ['t']), h[2] + h[1] + 1)

Sym = SymmetricFunctions(QQ['t'])
h = Sym.h()
a(dual_k_theoretic_h(2, 2, base_ring=QQ['t']), h[2] + 2*h[1] + 3)

Sym = SymmetricFunctions(QQ['t'])
h = Sym.h()
a(dual_k_theoretic_h(2, 3, base_ring=QQ['t']), h[2] + 3*h[1] + 6)

# h_[2,1](x, [1, 1])
Sym = SymmetricFunctions(QQ)
h = Sym.h()
a(dual_k_theoretic_h([2, 1], [1, 1]), h[1]**2 + h[1]*h[2] + 2*h[1] + h[2] + 1)

# h_[1, 2](x, [2, 3])
Sym = SymmetricFunctions(QQ)
h = Sym.h()
a(dual_k_theoretic_h([1, 2], [2, 3]), 3*h[1]**2 + h[1]*h[2] + 2*h[2] + 12*h[1] + 12)


# test dual grothendieck function
Sym = SymmetricFunctions(QQ)
h = Sym.h()
a(dual_grothendieck_function([2, 1]), h[1]*h[2] + h[2] - h[3])




# ALL DONE!
print('Testing completed successfully!', end='')
end_time = time.time()
elapsed_time = end_time - start_time
print(' Elapsed time = {}'.format(elapsed_time))
