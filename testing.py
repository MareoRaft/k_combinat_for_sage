# Tools for testing purposes.
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

