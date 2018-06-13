===========================
k_combinat_for_sage
===========================

k-Schur and other "k-math-objects" for combinatorics.  This code is meant to be used with SageMath, and will eventually be migrated into SageMath itself.  The functions in the code come from research by `Jennifer Morse <http://math.virginia.edu/people/jlm6cj/>`_ and her colleagues.  You will see in the documentation references to various research papers, which are the relevant topics.


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
  * k-irreducible-k-shapes
  * root ideals

etc.  For a full list of functions, read the documentation.


Install or upgrade
--------------------
::

	$ sage -pip install --upgrade k_combinat_for_sage


Usage
---------------
Put your code in a file such as `myscript.py`.  Import all the functions with::

	from k_combinat_for_sage.all import *

and then use functions at will, such as::

	sp = SkewPartition([[4,3,3,1], [1]])
	if is_symmetric(sp):
		print "it's symmetric, boogie woogie woogie!"

Alternatively, you can import each module you need under a name of your choice::

	from k_combinat_for_sage import partition as P
	from k_combinat_for_sage import skew_partition as SP
	from k_combinat_for_sage import k_shape as kS
	from k_combinat_for_sage import root_ideal as RI

and then use functions appropriately.  For example::

	sp = SkewPartition([[4,3,3,1], [1]])
	if SP.is_symmetric(sp):
		print "it's symmetric, boogie woogie woogie!"

Finally, you can run your file with sage::

	$ sage myscript.py
	it's symmetric, boogie woogie woogie!

