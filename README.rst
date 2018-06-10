===========================
k_combinat_for_sage
===========================

k-Schur and other "k-math-objects" for combinatorics.  This code is meant to be used with SageMath, and will eventually be migrated into SageMath itself.  The functions in the code come from research by `Jennifer Morse <http://math.virginia.edu/people/jlm6cj/>`_ and her colleagues.  You will see in the documentation references to various research papers, which are the relevant topics.


Quicklinks
--------------

  * `Documentation <https://mareoraft.github.io/morse-code/>`_
  * `proof of work <https://github.com/MareoRaft/morse-code/blob/master/src/proof_of_work.py>`_


Install
---------------
::

	$ sage -pip install k_combinat_for_sage


Usage
---------------
Put your code in a file such as `myscript.py`.  Import all the functions with::

	from k_combinat_for_sage.all import *
	
and then use functions at will, such as::

	sp = SkewPartition([[4,3,3,1], [1]])
	if SP.is_symmetric(sp):
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


Overview
---------------

This code computes combinatorial things such as

  * k-boundaries
  * k-rim
  * k-shape partitions
  * skew-linked diagrams
  * k-irreducible partitions
  * k-irreducible-k-shapes
  * root ideals

etc.


build
---------------

To build the documentation, comment out all sage-specific imports (such as `from sage.all import *`), cd into `doc` folder, and run::

	make html

our config file `doc/source/conf.py` may try to be similar to the official sage config file `src/doc/common/conf.py`.  Documentation can contain mathjax using backticks, for example, `\\sum \\frac{h}{l}`.

To automatically generate for gh pages, try "https://gist.github.com/brantfaircloth/791759" or "https://github.com/sphinx-doc/sphinx/issues/3382" or "https://daler.github.io/sphinxdoc-test/includeme.html"
