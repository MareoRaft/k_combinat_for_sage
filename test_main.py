#!/usr/bin/env sage
# A place to test my functions

from sage.all import *

from testing import *
from main import *
print('Sage loaded.  Testing...')


# test_right
sp = SkewPartition([[1], []])
c = SkewPartition_right(sp, 0)
a(c, 0)

sp = SkewPartition([[2], []])
c = SkewPartition_right(sp, 0)
a(c, 1)

sp = SkewPartition([[1, 1], []])
c = SkewPartition_right(sp, 0)
a(c, 0)

sp = SkewPartition([[2], [1]])
c = SkewPartition_right(sp, 0)
a(c, 1)

sp = SkewPartition([[2], [2]])
c = SkewPartition_right(sp, 0)
a(c, None)

sp = SkewPartition([[5, 5, 4, 2, 2],  [4, 3, 2]])
c = SkewPartition_right(sp, 2)
a(c, 3)


# test_left
sp = SkewPartition([[1], []])
c = SkewPartition_left(sp, 0)
a(c, 0)

sp = SkewPartition([[2], []])
c = SkewPartition_left(sp, 0)
a(c, 0)

sp = SkewPartition([[2], [1]])
c = SkewPartition_left(sp, 0)
a(c, 1)

sp = SkewPartition([[2], [2]])
c = SkewPartition_left(sp, 0)
a(c, None)

sp = SkewPartition([[5, 5, 4, 2, 2],  [4, 3, 2]])
c = SkewPartition_left(sp, 2)
a(c, 2)


# test_top
sp = SkewPartition([[1], []])
c = SkewPartition_top(sp, 0)
a(c, 0)

sp = SkewPartition([[2], []])
c = SkewPartition_top(sp, 0)
a(c, 0)

sp = SkewPartition([[1, 1], []])
c = SkewPartition_top(sp, 0)
a(c, 1)

sp = SkewPartition([[1, 1], [1]])
c = SkewPartition_top(sp, 0)
a(c, 1)

sp = SkewPartition([[2], [2]])
c = SkewPartition_top(sp, 0)
a(c, None)

sp = SkewPartition([[5, 5, 4, 2, 2],  [4, 3, 2]])
c = SkewPartition_top(sp, 2)
a(c, 2)

sp = SkewPartition([[5, 5, 4, 2, 2],  [4, 3, 2]])
c = SkewPartition_top(sp, 3)
a(c, 2)


# test_bottom
sp = SkewPartition([[1], []])
c = SkewPartition_bottom(sp, 0)
a(c, 0)

sp = SkewPartition([[2], []])
c = SkewPartition_bottom(sp, 0)
a(c, 0)

sp = SkewPartition([[1, 1], []])
c = SkewPartition_bottom(sp, 0)
a(c, 0)

sp = SkewPartition([[1, 1], [1]])
c = SkewPartition_bottom(sp, 0)
a(c, 1)

sp = SkewPartition([[2], [2]])
c = SkewPartition_bottom(sp, 0)
a(c, None)

sp = SkewPartition([[5, 5, 4, 2, 2],  [4, 3, 2]])
c = SkewPartition_bottom(sp, 2)
a(c, 2)

sp = SkewPartition([[5, 5, 4, 2, 2],  [4, 3, 2]])
c = SkewPartition_bottom(sp, 3)
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


# test_kShape_is_k_reducible_by_rectangle
s = Partition([1])
k = 1
(w,h) = (1,1)
a(kShape_is_k_reducible_by_rectangle(s, k, (w,h)), True)

s = Partition([2, 1])
k = 1
(w,h) = (1,1)
a(kShape_is_k_reducible_by_rectangle(s, k, (w,h)), True)

s = Partition([1, 1])
k = 2
(w,h) = (1,1)
a(kShape_is_k_reducible_by_rectangle(s, k, (w,h)), False)

s = Partition([1, 1])
k = 2
(w,h) = (1,2)
a(kShape_is_k_reducible_by_rectangle(s, k, (w,h)), True)

s = Partition([1, 1])
k = 2
(w,h) = (2,1)
a(kShape_is_k_reducible_by_rectangle(s, k, (w,h)), False)

s = Partition([2, 1, 1])
k = 2
(w,h) = (2,1)
a(kShape_is_k_reducible_by_rectangle(s, k, (w,h)), False)

s = Partition([2, 1, 1])
k = 2
(w,h) = (1,2)
a(kShape_is_k_reducible_by_rectangle(s, k, (w,h)), True)

s = Partition([2, 1, 1])
k = 2
(w,h) = (1,1)
a(kShape_is_k_reducible_by_rectangle(s, k, (w,h)), True)

s = Partition([2, 1, 1])
k = 3
(w,h) = (3,1)
a(kShape_is_k_reducible_by_rectangle(s, k, (w,h)), False)

s = Partition([2, 1, 1])
k = 3
(w,h) = (2,2)
a(kShape_is_k_reducible_by_rectangle(s, k, (w,h)), False)

s = Partition([2, 1, 1])
k = 3
(w,h) = (1,3)
a(kShape_is_k_reducible_by_rectangle(s, k, (w,h)), False)

s = Partition([2, 1, 1])
k = 3
(w,h) = (2,1)
a(kShape_is_k_reducible_by_rectangle(s, k, (w,h)), False)

s = Partition([2, 1, 1])
k = 3
(w,h) = (1,2)
a(kShape_is_k_reducible_by_rectangle(s, k, (w,h)), True)

s = Partition([3, 2, 1])
k = 3
(w,h) = (3,1)
a(kShape_is_k_reducible_by_rectangle(s, k, (w,h)), False)

s = Partition([3, 2, 1])
k = 3
(w,h) = (2,2)
a(kShape_is_k_reducible_by_rectangle(s, k, (w,h)), True)

s = Partition([3, 2, 1])
k = 3
(w,h) = (1,3)
a(kShape_is_k_reducible_by_rectangle(s, k, (w,h)), False)

s = Partition([3, 2, 1])
k = 3
(w,h) = (2,1)
a(kShape_is_k_reducible_by_rectangle(s, k, (w,h)), False)

s = Partition([3, 2, 1])
k = 3
(w,h) = (1,2)
a(kShape_is_k_reducible_by_rectangle(s, k, (w,h)), False)


# test_kShape_is_k_reducible2
s = Partition([1])
k = 1
a(kShape_is_k_reducible2(s, k), True)

s = Partition([2, 1])
k = 1
a(kShape_is_k_reducible2(s, k), True)

s = Partition([1, 1])
k = 2
a(kShape_is_k_reducible2(s, k), True)

s = Partition([2, 1, 1])
k = 2
a(kShape_is_k_reducible2(s, k), True)

s = Partition([2, 1, 1])
k = 3
a(kShape_is_k_reducible2(s, k), True)

s = Partition([3, 2, 1])
k = 3
a(kShape_is_k_reducible2(s, k), True)

s = Partition([5, 3, 2, 1, 1])
k = 4
a(kShape_is_k_reducible2(s, k), False)

s = Partition([5, 4, 2, 2, 1])
k = 4
a(kShape_is_k_reducible2(s, k), False)


# test_get_k_irreducible_k_shapes
ptns = get_k_irreducible_k_shapes(1)
a(ptns, [[]])

ptns = get_k_irreducible_k_shapes(2)
a(ptns, [[]])

ptns = get_k_irreducible_k_shapes(3)
a(ptns, [[], [1], [2, 1]])


# test_is_linked
# empty skew
sp = SkewPartition([[], []])
a(is_linked(sp), True)

# valid row shape but invalid col shape
sp = SkewPartition([[3, 2], [1, 0]])
a(is_linked(sp), False)



# ALL DONE!
print('Testing completed successfully!')
