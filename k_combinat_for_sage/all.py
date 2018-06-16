# -*- coding: utf-8 -*-
"""
This module contains all functionalities that are not already organized into the other files.  New functionalities written to the library often appear here, and eventually get organized into separate files.

REFERENCES:

.. [fun] `Raising operators and the Littlewood-Richardson polynomials <https://arxiv.org/pdf/1203.4729.pdf>`_.  Fun, Alex.
"""
from sage.all import *
# from sage.structure.unique_representation import UniqueRepresentation
# from sage.structure.parent import Parent
# from sage.categories.infinite_enumerated_sets import InfiniteEnumeratedSets

from partition import *
import partition as P
from skew_partition import *
import skew_partition as SP
from k_shape import *
import k_shape as kS
from root_ideal import *
import root_ideal as RI


# MAIN
def get_k_rectangles(k):
    # A __k-rectangle__ is a partition whose Ferrer's diagram is a rectangle whose largest hook-length is k.
    return [Partition([a] * b) for (a, b) in k_rectangle_dimension_list(k)]

def get_k_irreducible_partition_lists(k):
    # Returns: list of lists (instead of list of Partition objects)

    # Since there are n! such partitions, the big-O time can't be better than that.
    # We could have a yeild in the function to be an iterator.
    k = NonNegativeIntegerSemiring()(k)
    k_irr_ptns = [[]]
    # NO rows of length k
    for i in range(1, k):
        new_k_irr_ptns = []
        for ptn in k_irr_ptns:
            # at most i rows of length k-i where 1 <= i < k
            for num_rows in range(0, i+1):
                new_ptn = ptn + [k-i]*num_rows
                new_k_irr_ptns.append(new_ptn)
        k_irr_ptns = new_k_irr_ptns
    return k_irr_ptns

def get_k_irreducible_partitions(k):
    """ Given k, return the n! k-irreducible-partitions. """
    return [Partition(e) for e in get_k_irreducible_partition_lists(k)]

def size_to_num_linked_partition_self_pairs(n):
    # Given a natural number n, count how many partitions l of size n have the property that (l, l) has a corresponding linked-skew-diagram.
    ps = Partitions(n)
    count = 0
    for p in ps:
        try:
            row_col_to_skew_partition(p, p)
        except:
            pass
        else:
            count += 1
    return count

def print_sequence(func, num_terms=float('inf')):
    n = 0
    while n < num_terms:
        print('n={}\t{}=f(n)'.format(n, func(n)))

def size_to_k_shapes(n, k):
    """ Find all partitions of size n that are k-shapes. """
    return [ptn for ptn in Partitions(n) if is_k_shape(ptn, k)]

def size_to_num_k_shapes(n, k):
    return len(size_to_k_shapes(n, k))

def straighten(s, gamma):
    """ Perform Schur function straightening by the Schur straightening rule ([cat]_, Prop. 4.1).  Also known as the slinky rule.

    `s_\\gamma(\\mathbf{x}) = \\begin{cases}
        \\sgn(\\gamma+\rho) s_{\\text{sort}(\\gamma+\\rho) -\\rho}(\\mathbf{x}) & \\text{if $\\gamma + \\rho$ has distinct nonnegative parts,}\\
        0                                                          & \\text{otherwise,}
    \\end{cases}`

    where `\\rho=(\\ell-1,\\ell-2,\\dots,0)`, `\\text{sort}(\\beta)` denotes the weakly decreasing sequence obtained by sorting `\\beta`, and `\\sgn(\\beta)` denotes the sign of the (shortest possible) sorting permutation.

    EXAMPLES::

        sage: straighten([2, 1, 3])
        -s[2, 2, 2]
        # because s[2, 1, 3] := -s[2, 2, 2]
    """
    def has_nonnegative_parts(lis):
        return all(e >= 0 for e in lis)
    def has_distinct_parts(lis):
        return len(set(lis)) == len(lis)
    def number_of_noninversions(lis):
        num = 0
        for i in range(len(lis)):
            for j in range(i + 1, len(lis)):
                # i < j is already enforced
                if lis[i] < lis[j]:
                    num += 1
        return num
    rho = list(range(len(gamma) - 1, -1, -1))
    combined = [g + r for g, r in zip(gamma, rho)]
    if has_distinct_parts(combined) and has_nonnegative_parts(combined):
        sign = (-1)**number_of_noninversions(combined)
        sort_combined = reversed(sorted(combined))
        new_gamma = [sc - r for sc, r in zip(sort_combined, rho)]
        return sign * s(new_gamma)
    else:
        return 0


class ShiftingSequenceSpace():
    # A helper for ShiftingOperatorAlgebra
    def __init__(self, base=IntegerRing()):
        self.base = base
        # category = InfiniteEnumeratedSets()
        # Parent.__init__(self, category=category)

    def __contains__(self, seq):
        if not isinstance(seq, tuple):
            return False
        return not any(i not in self.base for i in seq)

    VALIDATION_ERROR_MESSAGE = 'Expected valid index (a tuple of {base}), but instead received {seq}.'
    def validate(self, seq):
        if not self.__contains__(seq):
            raise ValueError(self.VALIDATION_ERROR_MESSAGE.format(base=self.base, seq=seq))


class RaisingSequenceSpace(ShiftingSequenceSpace):
    # helper for RaisingOperatorAlgebra
    VALIDATION_ERROR_MESSAGE = 'Expected valid index (a tuple of {base} elements, where every partial sum is nonnegative and every total sum is 0), but instead received {seq}.'
    def __contains__(self, seq):
        # check that it is a shifting sequence
        if not ShiftingSequenceSpace.__contains__(self, seq):
            return False
        # check that every partial sum is nonnegative
        partial_sum = 0
        for term in seq:
            partial_sum += term
            if partial_sum < 0:
                return False
        # check that total sum is 0
        if partial_sum != 0:
            return False
        # finally, succeed
        return True


class ShiftingOperatorAlgebra(CombinatorialFreeModule):
    def __init__(self, base_ring=QQ['t'], prefix='S', basis_indecis=ShiftingSequenceSpace()):
        self._prefix = prefix
        self._base_ring = base_ring
        # a single basis index looks like (1, 0, -1, 2), for example
        self._basis_indecis = basis_indecis
        # category
        category = Algebras(self._base_ring.category()).WithBasis()
        category = category.or_subcategory(category)
        # init
        CombinatorialFreeModule.__init__(
            self,
            self._base_ring,
            self._basis_indecis,
            category=category,
            prefix=self._prefix,
            bracket=False)

    def __getitem__(self, seq):
        # seq should be a basis index
        self._basis_indecis.validate(seq)
        return self.basis()[seq]

    def _element_constructor_(self, seq):
        return self.__getitem__(seq)

    @cached_method
    def one_basis(self):
        # identity basis/index
        return tuple()

    def _repr_(self):
        return "Shifting Operator Algebra over {base_ring}".format(base_ring=self._base_ring)

    class Element(CombinatorialFreeModule.Element):
        """ element of a ShiftingOperatorAlgebra"""
        def indecis(self):
            return self.support()

        def index(self):
            if len(self) != 1:
                raise ValueError("This is only defined for basis elements.  For other elements, use indecis() instead.")
            return self.indecis()[0]

        def _mul_(self, other):
            def index_mul(index1, index2):
                max_len = max(len(index1), len(index2))
                index1 = index1 + (0,) * (max_len - len(index1))
                index2 = index2 + (0,) * (max_len - len(index2))
                return tuple(i1 + i2 for i1, i2 in zip(index1, index2))
            self_index_coeff_list = self.monomial_coefficients().items()
            other_index_coeff_list = other.monomial_coefficients().items()
            out_index_coeff_list = []
            for index1, coeff1 in self_index_coeff_list:
                for index2, coeff2 in other_index_coeff_list:
                    out_index_coeff = (index_mul(index1, index2), coeff1 * coeff2)
                    out_index_coeff_list.append(out_index_coeff)
            R = self.parent()
            out_list = [coeff * R(index) for index, coeff in out_index_coeff_list]
            out = R.sum(out_list)
            return out

        def __call__(self, operand):
            def raise_func(seq, operand):
                # seq is the index
                if isinstance(operand, (list, Composition, Partition)):
                    # pad sequence and operand with 0's
                    operand = operand + [0] * (len(seq) - len(operand))
                    seq = seq + (0,) * (len(operand) - len(seq))
                    # raise and drop
                    out = [v + s for v, s in zip(operand, seq)]
                else:
                    # it's some symmetric function basis element
                    parent_basis = operand.parent()
                    # process the vectors
                    dic = operand.monomial_coefficients()
                    assert len(dic) == 1
                    assert dic.values()[0] == 1
                    composition = dic.keys()[0]
                    out_composition = raise_func(seq, composition)
                    # TODO: check if out_composition is a valid partition.  when it's not, out becomes 0
                    # TODO: change below to 'out_partition' once validated
                    if parent_basis.__class__.__name__ == 'SymmetricFunctionAlgebra_schur_with_category':
                        out = straighten(parent_basis, out_composition)
                    else:
                        out = parent_basis(out_composition)
                return out
            def call_monomial(seq, coeff, operand, power=1):
                for _ in range(power):
                    operand = raise_func(seq, operand)
                return (operand, coeff)
            # start here
            if isinstance(operand, tuple):
                # if the operand is a tuple, perform __call__ on each piece (RECURSE)
                return tuple(self.__call__(op) for op in operand)
            # break into basis pieces
            index_coeff_list = self.monomial_coefficients().items()
            # perform raise_func on each piece
            out_list = [call_monomial(index, coeff, operand) for index, coeff in index_coeff_list]
            # recombine
            if isinstance(operand, list):
                out = out_list
            else:
                out = sum(coeff * mon for mon, coeff in out_list)
            return out

class RaisingOperatorAlgebra(ShiftingOperatorAlgebra):
    """
    We follow the following convention!:

    R((1, 0, -1)) is the raising operator that raises the first part by 1 and lowers the third part by 1.

    For a definition of raising operators, see [cat]_ Definition 2.1, but be wary that the notation is different there.  See :meth:`ij` for a way to create operators using the notation in the paper.

    If you do NOT want any restrictions on the allowed sequences, simply use 'ShiftingOperatorAlgebra' instead of 'RaisingOperatorAlgebra'.

    OPTIONAL ARGUMENTS:
    - ``base_ring`` -- (default ``QQ['t']``) the ring you will use on the raising operators.
    - ``prefix`` -- (default ``"R"``) the label for the raising operators.

    EXAMPLES::

        sage: R = RaisingOperatorAlgebra()
        sage: s = SymmetricFunctions(QQ['t']).s()
        sage: h = SymmetricFunctions(QQ['t']).h()

        sage: R((1, -1))
        R(1, -1)
        sage: R((1, -1))(s[5, 4])
        s[6, 3]
        sage: R((1, -1))(h[5, 4])
        h[6, 3]

        sage: (1 - R((1,-1))) * (1 - R((0,1,-1)))
        R() - R(0, 1, -1) - R(1, -1) + R(1, 0, -1)
        sage: ((1 - R((1,-1))) * (1 - R((0,1,-1))))(s[2, 2, 1])
        (-3*t-2)*s[] + s[2, 2, 1] - s[3, 1, 1] + s[3, 2]
    """
    def __init__(self, base_ring=QQ['t'], prefix='R'):
        ShiftingOperatorAlgebra.__init__(self,
            base_ring=base_ring,
            prefix=prefix,
            basis_indecis=RaisingSequenceSpace())

    def ij(self, i, j):
        """ Shorthand element constructor that allows you to create raising operators using the familiar `R_{ij}` notation found in [cat]_ Definition 2.1, with the exception that indecis here are 0-based, not 1-based.

        EXAMPLES::

            # create the raising operator which raises part 0 and lowers part 2 (indecis are 0-based)
            sage: R.ij(0, 2)
            R((1, 0, -1))

        """
        if not i in NonNegativeIntegerSemiring():
            raise ValueError('i must be a natural number.  You input i = {i}.'.format(i=i))
        if not j in NonNegativeIntegerSemiring():
            raise ValueError('j must be a natural number.  You input j = {j}.'.format(j=j))
        if not i < j:
            raise ValueError('Index j must be greater than index i.  You input (i, j) = ({i}, {j}).'.format(i=i, j=j))
        seq = [0] * (max(i, j) + 1)
        seq[i] = 1
        seq[j] = -1
        seq = tuple(seq)
        return self._element_constructor_(seq)


class HallLittlewoodVertexOperator:
    """
    Garsia's version of Jing's Hall-Littlewood vertex operators.  These are defined in equations 4.2 and 4.3 of [cat]_ and appear visually as a bold capital H.

    INPUTS:

    base_ring: (defaults to QQ['t']) the base ring to build the SymmetricFunctions upon.

    EXAMPLES::

        sage: H = HallLittlewoodVertexOperator
        sage: one = SymmetricFunctions(QQ['t']).hall_littlewood().Qp().one()
        sage: H([4, 1, 3])(one) == H(4)(H(1)(H(3)(one)))
        True

    """
    def __init__(self, composition, base_ring=QQ['t']):
        if composition in NonNegativeIntegerSemiring():
            self.composition = [composition]
        elif isinstance(composition, (list, Composition, Partition)):
            self.composition = composition
        else:
            raise ValueError('Bad composition.')
        self.base_ring = base_ring

    def __repr__(self):
        return '{}({})'.format(self.__class__.__name__, self.composition)

    def __call__(self, input_):
        gamma = self.composition
        # iterate
        for part in reversed(gamma):
            input_ = input_.hl_creation_operator([part])
        return input_


def compositional_hall_littlewood_Qp(gamma, base_ring=QQ['t']):
    """
    Given gamma, returns the compositional Hall-Littlewood polynomial `H_{\\gamma}(\\mathbf{x}; t)` in the Q' basis, as defined in [cat]_ section 4.4.

    If the composition gamma is a partition, this is just the Hall-Littlewood Q' polynomial.

    EXAMPLES::

        sage: hl = SymmetricFunctions(QQ['t']).hall_littlewood().Qp()
        sage: compositional_hall_littlewood_Qp([3, 3, 2]) == hl[3, 3, 2]
        True

    """
    sym = SymmetricFunctions(base_ring)
    hl = sym.hall_littlewood().Qp()
    H = HallLittlewoodVertexOperator
    return H(gamma)(hl.one())

def raising_root_ideal_operator(ri, t=1, base_ring=QQ['t']):
    """ Given a root ideal `ri = \\Phi` (and optionally a variable `t`), return the operator `\\Prod_{(i,j) \\in \\Phi} (1 - tR_{ij})`.
    """
    R = RaisingOperatorAlgebra(base_ring=base_ring)
    def prod(iterable):
        return reduce(R.Element._mul_, iterable, R.one())
    op = prod([1 - t*R.ij(ij) for ij in ri])
    return op

def indexed_root_ideal_to_catalan_function(ri, index, base_ring=QQ['t']):
    """
    INPUTS:

    ri: root_ideal
    index: composition that indexes the root ideal

    OUTPUTS:

    The catalan function
    """
    # setup
    hl = SymmetricFunctions(base_ring).hall_littlewood().Qp()
    t = base_ring.gen()
    # formula
    n = len(index)
    ri_complement = RI.complement(ri, n)
    op = raising_root_ideal_operator(ri_complement, t=t)
    cat_func = op(hl(index))
    return cat_func

def skew_partition_to_catalan_function(sp, base_ring=QQ['t']):
    """ Given a SkewPartition `sp = (\\lambda, \\mu)`, return the catalan function `H(\\Phi^+(sp); \\lambda)`.
    """
    ri = skew_partition_to_root_ideal(sp, type='max')
    rs = sp.row_lengths()
    return indexed_root_ideal_to_catalan_function(ri, rs, base_ring)

def row_and_column_lengths_to_catalan_function(row_lengths, column_lengths, base_ring=QQ['t']):
    """ Determine the skew partition `D` with row-shape `row_lengths` and column-shape `column_lengths`, and return the catalan function `H(\\Phi^+(D); row_lengths)`.
    """
    sp = SkewPartitions().from_row_and_column_length(row_lengths, column_lengths)
    return skew_partition_to_catalan_function(sp, base_ring)

def k_shape_to_catalan_function(p, k, base_ring=QQ['t']):
    """ Given `k` and a `k`-shape `p`, return the catalan function `H(\\Psi^+((rs(p),cs(p))), rs(p))`.
    """
    assert is_k_shape(p, k)
    rs = p.row_lengths()
    cs = p.column_lengths()
    return row_and_column_lengths_to_catalan_function(rs, cs, base_ring)

def k_plus_one_core_to_k_schur_function(p, k, base_ring=QQ['t']):
    # TODO: compare the performance of this function to existing k-schur function.
    assert is_k_core(p, k + 1)
    return k_shape_to_catalan_function(p, k, base_ring)

# def k_bdd_one_core_to_k_schur_function(p, k, base_ring=QQ['t']):
#     # TODO: compare the performance of this function to existing k-schur function.
#     assert is_k_core(p, k + 1)
#     return k_shape_to_catalan_function(p, k, base_ring)

class DoubleIntegers():
    # merely a helper for DoubleRing
    def __contains__(self, el):
        return el in IntegerRing() or el == '*'

# NEW IDEA.  Can we just use SYMMETRIC FUNCTIONS plugging in x = a and LOSING the SYMMETRIC PART?
class DoubleRing(CombinatorialFreeModule):
    # it's an algebra if you consider the a's and the decorative integers as 'separate'.
    # if you consider them all together, it's a ring
    def __init__(self, base_ring=IntegerRing(), prefix='a', basis_indecis=DoubleIntegers()):
        self._prefix = prefix
        self._base_ring = base_ring
        self._basis_indecis = basis_indecis
        # category
        # TODO: make commutative? make free?
        category = Algebras(self._base_ring.category()).WithBasis()
        category = category.or_subcategory(category)
        # init
        CombinatorialFreeModule.__init__(
            self,
            self._base_ring,
            self._basis_indecis,
            category=category,
            prefix=self._prefix,
            bracket=False)

    def __getitem__(self, index):
        assert index in self._basis_indecis
        return self.basis()[index]

    def _element_constructor_(self, index):
        return self.__getitem__(index)

    @cached_method
    def one_basis(self):
        # identity index
        return '*'

    def _repr_(self):
        return "DoubleRing over {base_ring}".format(base_ring=self._base_ring)

    class Element(CombinatorialFreeModule.Element):
        def indecis(self):
            return self.support()

        def index(self):
            if len(self) != 1:
                raise ValueError("This is only defined for basis elements.  For other elements, use indecis() instead.")
            return self.indecis()[0]

        # def _mul_(self, other):
            # TODO: make it free
            # we could change basis elements to a_dic where dic keys is an increasing list of integers and values are multiplicities
            # this would be free and have multiplication and be commutative.


