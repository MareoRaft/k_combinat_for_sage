# -*- coding: utf-8 -*-
"""
This module contains all functionalities that are not already organized into the other files.  New functionalities written to the library often appear here, and eventually get organized into separate files.
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
    k = NN(k)
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

def n_to_number_of_linked_partition_self_pairs(n):
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

def n_to_k_shapes(n, k):
    """ Given n, find all partitions of size n that are k-shapes. """
    return [ptn for ptn in Partitions(n) if is_k_shape(ptn, k)]

def n_to_num_k_shapes(n, k):
    return len(n_to_k_shapes(n, k))

def n_to_k_shape_boundaries(n, k):
    # Given n, find all k-boundaries of all k-shapes of size n.
    return [ptn.k_boundary(k) for ptn in Partitions(n) if is_k_shape(ptn, k)]

def n_to_symmetric_k_shape_boundaries(n, k):
    k_shape_boundaries = n_to_k_shape_boundaries(n, k)
    return [ks for ks in k_shape_boundaries if kS.is_symmetric(ks)]

def n_to_num_symmetric_k_shape_boundaries(n, k):
    return len(n_to_symmetric_k_shape_boundaries(n, k))


class ShiftingSequenceSpace():
    # A helper for ShiftingOperatorAlgebra
    def __init__(self, base=ZZ):
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
                if isinstance(operand, list) or isinstance(operand, Partition):
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

    EXAMPLE::

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

        EXAMPLE::

            # create the raising operator which raises part 0 and lowers part 2 (indecis are 0-based)
            sage: R.ij(0, 2)
            R((1, 0, -1))

        """
        seq = [0] * (max(i, j) + 1)
        seq[i] = 1
        seq[j] = -1
        seq = tuple(seq)
        return self._element_constructor_(seq)


def straighten(s, gamma):
    """ Perform Schur function straightening by the Schur straightening rule ([cat]_, Prop. 4.1).

    `s_\\gamma(\\mathbf{x}) = \\begin{cases}
        \\sgn(\\gamma+\rho) s_{\\text{sort}(\\gamma+\\rho) -\\rho}(\\mathbf{x}) & \\text{if $\\gamma + \\rho$ has distinct nonnegative parts,}\\
        0                                                          & \\text{otherwise,}
    \\end{cases}`

    where `\\rho=(\\ell-1,\\ell-2,\\dots,0)`, `\\text{sort}(\\beta)` denotes the weakly decreasing sequence obtained by sorting `\\beta`, and `\\sgn(\\beta)` denotes the sign of the (shortest possible) sorting permutation.

    EXAMPLE::

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



class HallLittlewoodVertexOperator:
    """
    Garsia's version of Jing's Hall-Littlewood vertex operators.  These are defined in equations 4.2 and 4.3 of [cat]_ and appear visually as a bold capital H.

    INPUTS:

    base_ring: (defaults to QQ['t']) the base ring to build the SymmetricFunctions upon.

    EXAMPLE::

        sage: H = HallLittlewoodVertexOperator
        sage: one = SymmetricFunctions(QQ['t']).hall_littlewood().Qp().one()
        sage: H([4, 1, 3])(one) == H(4)(H(1)(H(3)(one)))
        True

    """
    def __init__(self, composition, base_ring=QQ['t']):
        if composition in NN:
            self.composition = [composition]
        elif isinstance(composition, list) or isinstance(composition, Partition):
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

    EXAMPLE::

        sage: hl = SymmetricFunctions(QQ['t']).hall_littlewood().Qp()
        sage: compositional_hall_littlewood_Qp([3, 3, 2]) == hl[3, 3, 2]
        True

    """
    sym = SymmetricFunctions(base_ring)
    hl = sym.hall_littlewood().Qp()
    H = HallLittlewoodVertexOperator
    return H(gamma)(hl.one())

def indexed_root_ideal_to_catalan_function(ri, index, base_ring=QQ['t']):
    """
    INPUTS:

    ri: root_ideal
    index: composition that indexes the root ideal

    OUTPUTS:

    The catalan function
    """
    # setup
    sym = SymmetricFunctions(base_ring)
    hl = sym.hall_littlewood().Qp()
    R = RaisingOperatorAlgebra()
    t = base_ring.gen()
    def prod(iterable):
        return reduce(R.Element._mul_, iterable, R.one())
    # formula
    n = len(index)
    ri_complement = RI.complement(ri, n)
    op = prod([1 - t*R.ij(ij) for ij in ri_complement])
    cat_func = op(hl(index))
    return cat_func

