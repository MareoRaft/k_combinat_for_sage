# MorseCode for SageMath

Code for Jennifer Morse.

I am working with [Jennifer Morse](http://math.virginia.edu/people/jlm6cj/) over the summer (and probably for a few years).  I am writing code with the goal of getting it added to the Sagemath project.  This code computes combinatorial things such as

  * k-boundaries
  * k-rim
  * k-shape partitions
  * skew-linked diagrams
  * k-irreducible partitions
  * k-irreducible-k-shapes
  * root ideals

etc.

## build
To build the documentation, comment out all sage-specific imports (such as `from sage.all import *`), cd into `doc` folder, and run:

    make html

our config file `doc/source/conf.py` may try to be similar to the official sage config file `src/doc/common/conf.py`.  Documentation can contain mathjax using backticks, for example, `\\sum \\frac{h}{l}`.

To automatically generate for gh pages, try "https://gist.github.com/brantfaircloth/791759" or "https://github.com/sphinx-doc/sphinx/issues/3382" or "https://daler.github.io/sphinxdoc-test/includeme.html"
