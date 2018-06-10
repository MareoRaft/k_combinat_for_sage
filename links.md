Notes/links just for me

https://ask.sagemath.org/question/9302/import-sage-packages-in-python/

see
https://ask.sagemath.org/question/33954/can-i-create-a-sage-file-and-import-it-as-a-python-module/
for explanation of import_statements and preparse

classcall private explained:
https://doc.sagemath.org/html/en/reference/structure/sage/structure/unique_representation.html

sage track ticket:
https://trac.sagemath.org/ticket/25295
https://trac.sagemath.org/ticket/25298
https://trac.sagemath.org/ticket/25538


sage combinat things to take a look out and see what's there:
ncsf_qsym
sf/k_dual
sf/macdonald
sf/new_kschur
k_tableau



things that definitely exist in sage
------------------------------------
  * dual k-Schur basis
  * strong and weak k-tableau
  * SkewTableau
  * SkewPartition (ferrers diagram of the skew shape)
  * row-shape and col-shape exist in SkewPartition as row_lengths and column_lengths
  * "cells": We use "cells" as coordinates.  A "cell" is an ordered pair (i, j).  No actual class.
  * hook_length and hook_lengths methods in Partition
  * leg_length and arm_length
  * k_boundary and k_interior already exist in Partition!
  * cell_poset exists, which is just like root_ideal way of listing cells in a partition, except it is left-justified instead of right-justified.  IN FACT, it allows you to set the 'orientation' as NE (northeast), so we may get rid of some of our root_ideal code.
  * corners.  given a ptn, then ptn.corners() is the corners of the Young diagram!
  * k_conjugate (whatever that is)
  * k_atom
  * k_split


things mentioned in sage (and might exist)
------------------------------------------
  * k-bounded partitions
  * kshapes mentioned in k_tableau and skew_tableau



things that definitely do not exist in sage
-------------------------------------------
  * 



plan
----
If I need to see whether an object (for example a skew shape) qualifies to be an object of a child class (for example a k-boundary), there are three options:

(1) Create `is_k_boundary` as a method of the skew shape class.
(2) Create a separate `is_k_boundary` function in my module for k-things.  The point being that we don't need to bloat the skew shape class.  (Of course, kBoundary class will use is_k_boundary function itself).
(3) Just do the checking in the `__init__` method of the kBoundary class ONLY.  If the user wants to see if a skew shape is a kBoundary, they must try to create a new kBoundary(skew_shape) object and catch the error.

Planning to go with (1) because it seems consistent with what I'm seeing: In the Tableau class, I see methods like is_standard, is_k_tableau, and is_key_tableau, two of which I know have their own class.






Questions
-----------
I think a lot of these objects will be SkewPartitions or ferrers diagrams, but then should have a .to_hook_length_SkewTableau method which gives you back the shape with hook lengths filled in.


build
---------------

To build the documentation, comment out all sage-specific imports (such as `from sage.all import *`), cd into `doc` folder, and run::

	make html

our config file `doc/source/conf.py` may try to be similar to the official sage config file `src/doc/common/conf.py`.  Documentation can contain mathjax using backticks, for example, `\\sum \\frac{h}{l}`.

To automatically generate for gh pages, try "https://gist.github.com/brantfaircloth/791759" or "https://github.com/sphinx-doc/sphinx/issues/3382" or "https://daler.github.io/sphinxdoc-test/includeme.html"


listing in PyPI
------------------
Trying to add repo to PyPI.  First tried
https://peterdowns.com/posts/first-time-with-pypi.html
, now trying
https://stackoverflow.com/questions/45207128/failed-to-upload-packages-to-pypi-410-gone/45209514#45209514






python setup.py sdist
twine upload dist/*
