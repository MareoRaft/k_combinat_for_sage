from sage.all import *
from sage.combinat.skew_partition import SkewPartition
# from sage.combinat.

# HELPERS (only exist as helper functions for other things):
def is_weakly_decreasing(li):
	return all(li[i] >= li[i+1] for i in range(len(li)-1))


# MAIN:


class SkewPartition2 (SkewPartition):
	def is_skew_linked_diagram(self):
		"""
		A __skew-linked diagram__ is a skew-shape `s` where both the row-shape and column-shape of `s` are partitions.
		"""
		return is_weakly_decreasing(self.row_lengths()) and is_weakly_decreasing(self.column_lengths())

	# DOESN"T MAKE SENSE BECAUSE OUR SKEWPARTITION ARE JUST OUTER AND INNER PART ALREADY DEFINED.
	# def minimal_containing_partition(self):
	# 	""" Returns the smallest possible partition that could contain this skew shape.
	# 	Simply go down and left of all cells in skew_shape
	# 	"""
	# 	current_row_len = 0
	# 	partition_reversed = []
	# 	for row in reversed(skew_shape.rows()):
	# 		current_row_len = max(current_row_len, rightmost(row))
	# 		partition_reversed.append(current_row_len)
	# 	partition = list(reversed(partition_reversed))
	# 	return partition


class kBoundary (SkewPartition2):
	"""
	Given a partition \lambda and a positive integer k, the __k-boundary__ of \lambda is the skew-shape obtained from the shape of \lambda by removing all cells of hook-length greater than k.
	"""
	def __init__(self, l, k):
		# NOTE: THIS FUNCTION IS REDUNDANT WITH Permutation.k_boundary AND THIS SHOULD BE ADDRESSED
		"""
		l: the partition
		k: the largest allowed hook length
		"""
		# to make a SkewPartition, we must calculate the inner and outer partitions
		outer_partition = l.to_list()
		# we could make more efficient, but it's simple to get the tableau of all hook lengths and pick out the inner partition from it
		inner_partition = []
		for row_hook_lengths in l.hook_lengths():
			inner_row_length = len([x for x in row_hook_lengths if x > k])
			inner_partition.append(inner_row_length)
		# finally, make the SkewPartition
		SkewPartition2.__init__(self, [outer_partition, inner_partition])

	def partition(self):
		""" Return the partition whose k-boundary is self. """



"""
Given a skew-linked diagram, is it a k-boundary?  (That is, does there exist some partition which - when cells of hook-length > k are removed - becomes the skew-linked diagram.)
"""

def k_boundary_to_partition(skew_shape, strict=True):
	"""
	skew_shape: The skew-shape (a k-boundary) to find the corresponding partition for.
	strict: If true, assert that the skew-shape really is a k-boundary.
	"""
	if strict:
		assert is_k_boundary(skew_shape)
	# return minimal_containing_partition()
	return skew_shape.outer()

def is_k_boundary(skew_shape):
	if k == 0:
		return is_empty_tableau(skew_shape)
	else:
		"""We go down and left of each cell to create the only possible partition that could have led to this potential k-boundary

		(Any other partition containing this skew_shape would necessarily have a northeast corner that the skew_shape does *not have*.  But in order for the skew-shape to be a k-boundary, it *must have* that northeast corner.)
		"""
		l = k_boundary_to_partition(skew_shape, strict=False)
		"""now that we have the partition, we simply compute it's hook-length for each cell and verify that for each cell of values k or less, it appears in the skew_shape"""
		for cell in l:
			if hook_length(l, cell) <= k:
				if cell not in skew_shape:
					return False
			else:
				if cell in skew_shape:
					return False
		return True

"""
Definition: A partition is __k-irreducible__ if its shape has at most k-i rows of length i for all 1 \leq i < k, and no rows of length \geq k.
"""

def is_k_irreducible(l, k):
	"""
	l: a partition
	k: the 'k_irreducible' k
	"""
	""" if there's more than k-i rows of length i, it fails the condition """
	assert is_partition(l) # we assume l is weakly decreasing
	if k < 0:
		raise ValueError('Nonsensical choice for k.')
	if k == 0:
		# do anything speciall???
		# return is_empty(l) ?

	l_reversed = list(reversed(l))
	i = 1
	num_rows_of_len_i = 0
	index = 0
	while True:
		if index == len(l_reversed) or l_reversed[index] > i:
			# check for failure, move on to next i
			if num_rows_of_len_i > k-i:
				return False
			i += 1
			num_rows_of_len_i = 0
		else:
			# add to the count, keep going
			assert l_reversed[index] == i # sanity check
			num_rows_of_len_i += 1
			index += 1
	return True

def get_k_irreducible_partitions(k):
	""" Given k, output a list of all k-irreducible partitions """





