Notes/links just for me


classcall private explained:
https://doc.sagemath.org/html/en/reference/structure/sage/structure/unique_representation.html

sage trac ticket:
k_stuff:
https://trac.sagemath.org/ticket/25295
question:
https://trac.sagemath.org/ticket/25298
hall inner product:
https://trac.sagemath.org/ticket/25538
h/s/hl_creation_operator input behavior:
https://trac.sagemath.org/ticket/25734

ask sage stuff:
where is sage and sage's python installed:
http://ask.sagemath.org/question/34337/where-is-sage-installed/
codomain could not be determined error:
https://ask.sagemath.org/question/42732/codomain-could-not-be-determined/
importing sage into python:
https://ask.sagemath.org/question/9302/import-sage-packages-in-python/
explanation of import_statements and preparse:
https://ask.sagemath.org/question/33954/can-i-create-a-sage-file-and-import-it-as-a-python-module/

how to build sage:
possibly need to cd into sage repo
sage -b
OR
MAKE='make -jNUM' make
OR
make build
where NUM is number of threads you wan't to devote to installation
see
https://doc.sagemath.org/html/en/installation/source.html#step-by-step-installation-procedure
for more

ICERM:
https://app.icerm.brown.edu/Cube/#!status



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
  * hl_creation_operator(nu, t=None) -- supposedly more general than Jing's Hall Littlewood operators!  Reading https://core.ac.uk/down"If k=1, this is Garsia's version of JHLVO".  These use the Q' basis of hall littlewood.
  * strong tableau: https://doc.sagemath.org/html/en/reference/combinat/sage/combinat/k_tableau.html
  * I think an SMT is just a StandardTableau with the standard weight 1111111
  * a Core object is a k-core and located at sage.combinat.core.Core



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
To build, make sure you have committed things, saved files, and then run

  ./docs/build_docs.py

To automatically generate for gh pages, try "https://gist.github.com/brantfaircloth/791759" or "https://github.com/sphinx-doc/sphinx/issues/3382" or "https://daler.github.io/sphinxdoc-test/includeme.html"


listing in PyPI
------------------
Trying to add repo to PyPI.  First tried
https://peterdowns.com/posts/first-time-with-pypi.html
, now trying
https://stackoverflow.com/questions/45207128/failed-to-upload-packages-to-pypi-410-gone/45209514#45209514




https://doc.sagemath.org/html/en/thematic_tutorials/tutorial-implementing-algebraic-structures.html
read also:
https://doc.sagemath.org/html/en/thematic_tutorials/coercion_and_categories.html#coercion-and-categories
and now
https://doc.sagemath.org/html/en/reference/categories/sage/categories/algebras_with_basis.html#sage.categories.algebras_with_basis.AlgebrasWithBasis
to get the generators:
A.algebra_generators()
a free algebra:
sage.categories.examples.algebras_with_basis.FreeAlgebra_with_category
sage.categories.examples.algebras_with_basis.FreeAlgebra(R, alphabet=('a', 'b', 'c'))

magmatic:
C = sage.categories.magmatic_algebras.MagmaticAlgebras.WithBasis(QQ)

there is also a plain ole 'FreeAlgebra' category in sage.all
See :mod:`~sage.algebras.free_algebra` for examples and corner cases.


we need to figure out how to make a direct sum in sage.  It may be overfill, but we could always use CombinatorialFreeModule again, this time in the category of groups or modules, with basis [0, 1, 2, 3, 4, ...] and over the ring ZZ.

That would give just what we want.


A workaround for sphinx import errors (NEVERMIND.  IT STILL TRIES TO IMPORT SAGE).  Create a file called sage_imports.  Import all sage things in THAT file.  Other files simply
from sage_imports import *
This would fool sphinx since sage_imports is in the project. NEVERMIND.  IT STILL TRIES TO IMPORT SAGE

new idea, have a fake imports file that defines each of the variables:
"
class UniqueRepresentation:
  pass
class ZZ:
  pass
"
and see if it fools sphinx.





documentation:
look at
https://groups.google.com/forum/#!topic/sage-support/oJhV0fLdv7s
for help


SSL error:
------------
1. try injecting SSL certificates following:
https://stackoverflow.com/questions/19377045/pip-cert-failed-but-curl-works#19398611
try the --trusted-host option in the pip install line

2. try updating openssl:
/usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"
brew upgrade openssl
and brew may tell you to then add something like
export LDFLAGS="-L/usr/local/opt/openssl@1.1/lib"
export CPPFLAGS="-I/usr/local/opt/openssl@1.1/include"
export PKG_CONFIG_PATH="/usr/local/opt/openssl@1.1/lib/pkgconfig"
or
echo 'export PATH="/usr/local/opt/openssl/bin:$PATH"' >> ~/.bash_profile
export LDFLAGS="-L/usr/local/opt/openssl/lib"
export CPPFLAGS="-I/usr/local/opt/openssl/include"
as according to
https://stackoverflow.com/questions/41328451/ssl-module-in-python-is-not-available-when-installing-package-with-pip3#42328420
try the --trusted-host option in the pip install line

3. try downloading CORRECT version of sage for the OS.  We need the OS X El Capitan (10.11.6) version of sage
trying sage 8.2 tar.bz2
trying sage 7.6 app

4. This is exactly the issue in ask sage:
https://ask.sagemath.org/question/38746/sage-pip-not-compatible-with-pypi/
suggested solution:
sage -i openssl
sage -f python2
git clone https://github.com/pyca/pyopenssl.git
cd pyopenssl
sage -pip install -e





todo:
change example from is_symmetric to SOMETHING ELSE


SageDays things to cover, sage days, questions
---------------------------------------------------
  1. `sage -pip install` is *not* working by default due to a common SSL issue.  This should somehow be patched with a more robust sage installation.
  2. The 'codomain' error issue above.
  3. The classcall thing for super rsk (attempt this on your own first)
  4. Roadmap for migrating k_combinat_for_sage into sage
  5. Documentation searching needs to be more effective
  6. subclass issues.  Why not Core < Partition < Composition ?

