#!/usr/bin/env sage
from sage.all import *


def get_k_irreducible_partition_lists(k):
	"""matt
	Since there are n! such partitions, the big-O time can't be better than that.
	We could have a yeild in the function to be an iterator.
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
