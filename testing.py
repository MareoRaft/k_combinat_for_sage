# Tools for testing purposes.
def a(one, two):
    """ Replaces "assert" statements.  Tells you what the input what is the assertion fails. """
    if one == two:
        pass
    else:
        msg = "\nExpected: {}\nGot: {}".format(two, one)
        raise AssertionError(msg)
