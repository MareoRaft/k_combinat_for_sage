r"""
A place for things that may be useful in core sage (not specific to k-combinatorics)

AUTHORS:

- Matthew Lancellotti (2018): Initial version
"""

#*****************************************************************************
#  Copyright (C) 2018 Matthew Lancellotti <mvlancellotti@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.all import *
# ^*^ sphinx insert ^*^


# class InfiniteDimensionalFreeRing (CommutativeRing, InfiniteDimensionalFreeAlgebra):
# 	pass

# base_ring=IntegerRing()
# algebras = Algebras(base_ring.category()).WithBasis()
# commutative_rings = CommutativeRings()
# F = ForgetfulFunctor(algebras, commutative_rings)
# InfiniteDimensionalFreeRing = F(InfiniteDimensionalFreeAlgebra)

# idea: coerce

# idea: manually change the category
# found '_init_category_' '_initial_coerce_list' '_initial_convert_list' '_unset_category' 'category' 'categories' 'coerce' 'hom' 'is_ring'

# Let us declare a coercion from `\ZZ[x]` to `\ZZ[z]`::
#  |
#  |                  sage: Z.<z> = ZZ[]
#  |                  sage: phi = Hom(X, Z)(z)
#  |                  sage: phi(x^2+1)
#  |                  z^2 + 1
#  |                  sage: phi.register_as_coercion()
#  |
#  |              Now we can add elements from `\ZZ[x]` and `\ZZ[z]`, because
#  |              the elements of the former are allowed to be implicitly
#  |              coerced into the later::
#  |
#  |                  sage: x^2 + z
#  |                  z^2 + z

# idea: patch SymmetricFunctions to accept Algebras, not just 'commutative rings'
