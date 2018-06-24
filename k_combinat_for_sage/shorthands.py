""" For the convenience of a user using this library heavily, I have provided a bunch of truncated function names.  If you wish to use these 'shorthand' function names, simply::

	from k_combinat_for_sage.shorthands import *

"""
from all import *
# ^*^ sphinx insert ^*^

# pre-initialized useful variables
Sym = SymmetricFunctions(QQ['t'])
Sym.inject_shorthands() # creates s, h, etc

# shorter names for certain functions
dual_k_theoretic_h = dual_k_theoretic_homogeneous
Kh = dual_k_theoretic_h
chl = compositional_hall_littlewood_Qp
H = HallLittlewoodVertexOperator
R = RaisingOperatorAlgebra
st = straighten
double_h = double_homogeneous
double_s = double_schur
