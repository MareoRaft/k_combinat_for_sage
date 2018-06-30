# Tools for testing purposes.
from sage.all import *

def is_False_or_None(x):
	return x is None or x == False

def a(one, two=None, func=None):
    """ Replaces "assert" statements.  Tells you what the input what is the assertion fails. """
    if func is not None:
    	if func(one):
    		pass
    	else:
	        msg = "\nExpected: {}\nGot: {}".format(func, one)
	        raise AssertionError(msg)
    else:
	    if one == two:
	        pass
	    else:
	        msg = "\nExpected: {}\nGot: {}".format(two, one)
	        raise AssertionError(msg)

def partition_tuple(*args):
	# usage partition_tuple([2, 1], [3]) == (Partition([2, 1]), Partition([3]))
	return tuple([Partition(p) for p in args])
