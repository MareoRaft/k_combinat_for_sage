===========================
k_combinat_for_sage
===========================

k-Schur and other "k-math-objects" for combinatorics.  This code is meant to be used with ![SageMath](https://www.sagemath.org/), and parts of it have already been migrated into SageMath itself.  The functions in the code come from research by `Jennifer Morse <http://math.virginia.edu/people/jlm6cj/>`_ and her colleagues.  You will see in the documentation references to various research papers, which are the relevant topics.


Quicklinks
--------------

  * `documentation <https://mareoraft.github.io/k_combinat_for_sage/>`_
  * `proof of work <https://github.com/MareoRaft/k_combinat_for_sage/blob/master/k_combinat_for_sage/proof_of_work.py>`_


Contents
---------------

This code computes combinatorial things such as

  * k-boundaries
  * k-rim
  * k-shape partitions
  * skew-linked diagrams
  * k-irreducible partitions
  * irreducible k-shapes
  * root ideals
  * raising root operators
  * catalan functions
  * dual symmetric functions
  * double symmetric functions

etc.  For a full list of functions, read the documentation.


Install or upgrade
--------------------
::

	$ sage -pip install --upgrade k_combinat_for_sage

If you get an SSL error, you can download the latest version manually at `PyPI <https://pypi.org/project/k-combinat-for-sage/#files>`_ and then install/upgrade with :code:`sage -pip install --upgrade /path/to/downloaded/file`.  If you want to fix the SSL issue permanently, look `here <https://ask.sagemath.org/question/38746/sage-pip-not-compatible-with-pypi/>`_.


Usage
---------------
Put your code in a file such as `myscript.py`.  Import all the functions with::

	from k_combinat_for_sage.all import *

and then use functions at will, such as::

	p = Partition([4, 3, 3, 1])
	if is_symmetric(p):
		print("it's symmetric, boogie woogie woogie!")

Finally, you can run your file with sage::

	$ sage myscript.py
	it's symmetric, boogie woogie woogie!

