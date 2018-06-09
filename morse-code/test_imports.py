class k_boundary ():
	def __init__(l, k):
		# l: the partition
		# k: the largest allowed hook length
		self.super(partition=l) # effectively s = l.shape()
		s = self.shape
		for cell in s.cells():
			if s.hook_length(cell) > k:
				s.remove(cell)

class Fun ():
	def __init__(a, b):
		print(a + b)


def fun(a, b):
	return a + b
