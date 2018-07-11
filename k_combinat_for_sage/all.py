# -*- coding: utf-8 -*-
r"""
This module contains all functionalities that are not already organized into the other files.  New functionalities written to the library often appear here, and eventually get organized into separate files.

REFERENCES:

.. [Fun] `Raising operators and the Littlewood-Richardson polynomials <https://arxiv.org/pdf/1203.4729.pdf>`_.  Fun, Alex.
.. [LN] `Finite sum Cauchy identity for dual Grothendieck polynomials <https://projecteuclid.org/download/pdf_1/euclid.pja/1407415930>`_.

"""
from sage.all import *

from core import *
from partition import *
import partition as P
from skew_partition import *
import skew_partition as SP
from k_shape import *
import k_shape as kS
from root_ideal import *
import root_ideal as RI
from strong_marked_tableau import *
# ^*^ sphinx insert ^*^


# HELPERS


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
    r""" Given k, return the n! k-irreducible-partitions. """
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
    r""" Find all partitions of size n that are k-shapes. """
    return [ptn for ptn in Partitions(n) if is_k_shape(ptn, k)]

def size_to_num_k_shapes(n, k):
    return len(size_to_k_shapes(n, k))

def straighten(s, gamma):
    r""" Perform Schur function straightening by the Schur straightening rule ([cat]_, Prop. 4.1).  Also known as the slinky rule.

    .. MATH::

        s_{\gamma}(\mathbf{x}) = \begin{cases}
            \text{sgn}(\gamma+\rho) s_{\text{sort}(\gamma+\rho) -\rho}(\mathbf{x}) & \text{if } \gamma + \rho \text{ has distinct nonnegative parts,} \\
            0 & \text{otherwise,}
        \end{cases}

    where `\rho=(\ell-1,\ell-2,\dots,0)`, `\text{sort}(\beta)` denotes the weakly decreasing sequence obtained by sorting `\beta`, and `\text{sgn}(\beta)` denotes the sign of the (shortest possible) sorting permutation.

    EXAMPLES::

        sage: s = SymmetricFunctions(QQ).s()
        sage: straighten(s, [2, 1, 3])
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
    if s.__class__.__name__ in ('SymmetricFunctionAlgebra_monomial_with_category', 'SymmetricFunctionAlgebra_dual_with_category'):
        raise NotImplemented('Straightening does not exist (that i know of) for the monomial basis or the forgotten/dual basis.')
    elif s.__class__.__name__ == 'HallLittlewood_qp_with_category':
        # there is no straightening that i know of, except that for raising operators, invalid partitions evaluate to 0
        if has_nonnegative_parts(gamma) and is_weakly_decreasing(gamma):
            return s(gamma)
        else:
            return 0
    elif s.__class__.__name__ in ('SymmetricFunctionAlgebra_homogeneous_with_category', 'SymmetricFunctionAlgebra_elementary_with_category', 'SymmetricFunctionAlgebra_power_with_category', 'SymmetricFunctionAlgebra_witt_with_category'):
        new_gamma = list(reversed(sorted(gamma)))
        if has_nonnegative_parts(new_gamma):
            return s(new_gamma)
        else:
            return 0
    elif s.__class__.__name__ == 'SymmetricFunctionAlgebra_schur_with_category':
        rho = list(range(len(gamma) - 1, -1, -1))
        combined = [g + r for g, r in zip(gamma, rho)]
        if has_distinct_parts(combined) and has_nonnegative_parts(combined):
            sign = (-1)**number_of_noninversions(combined)
            sort_combined = reversed(sorted(combined))
            new_gamma = [sc - r for sc, r in zip(sort_combined, rho)]
            return sign * s(new_gamma)
        else:
            return 0
    else:
        raise ValueError("The input parameter 's' should be a symmetric function basis.  For example, 's = SymmetricFunctions(QQ).s(); straighten(s, [2, 1, 3])', or one could use 'h' instead of 's'.")


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
    def __init__(self, base_ring=QQ['t'], prefix='S', basis_indices=ShiftingSequenceSpace()):
        self._prefix = prefix
        self._base_ring = base_ring
        # a single basis index looks like (1, 0, -1, 2), for example
        self._basis_indices = basis_indices
        # category
        category = Algebras(self._base_ring.category()).WithBasis()
        category = category.or_subcategory(category)
        # init
        CombinatorialFreeModule.__init__(
            self,
            self._base_ring,
            self._basis_indices,
            category=category,
            prefix=self._prefix,
            bracket=False)

    def __getitem__(self, seq):
        # seq should be a basis index
        self._basis_indices.validate(seq)
        return self.basis()[seq]

    def _element_constructor_(self, seq):
        return self.__getitem__(seq)

    @cached_method
    def one_basis(self):
        # identity basis/index
        return tuple()

    def _repr_(self):
        return "Shifting Operator Algebra over {base_ring}".format(base_ring=self._base_ring)

    def product_on_basis(self, index1, index2):
        # pad with 0's
        max_len = max(len(index1), len(index2))
        index1 = index1 + (0,) * (max_len - len(index1))
        index2 = index2 + (0,) * (max_len - len(index2))
        # add the vectors
        index_product = tuple(i1 + i2 for i1, i2 in zip(index1, index2))
        return self.__getitem__(index_product)

    class Element(CombinatorialFreeModule.Element):
        r""" element of a ShiftingOperatorAlgebra"""
        def indices(self):
            return self.support()

        def index(self):
            if len(self) != 1:
                raise ValueError("This is only defined for basis elements.  For other elements, use indices() instead.")
            return self.indices()[0]

        def _call_basis_on_index(self, seq, index):
            # a 'seq' is the sequence of the BASIS element to act WITH.
            # an 'index' is a sequence (typically a composition or a partition) that we act UPON.
            assert is_sequence(index)
            # pad sequence and index with 0's
            index = index + [0] * (len(seq) - len(index))
            seq = seq + (0,) * (len(index) - len(seq))
            # raise and drop
            return [v + s for v, s in zip(index, seq)]

        def __call__(self, operand):
            def raise_func(seq, operand):
                # seq is the index
                if is_sequence(operand):
                    return self._call_basis_on_index(seq, operand)
                else:
                    # it's some symmetric function basis element
                    parent_basis = operand.parent()
                    # process the vectors
                    dic = operand.monomial_coefficients()
                    assert len(dic) == 1
                    (composition, coeff) = dic.items()[0] # occasionally a coefficient can show up (not cool, so consider the inclusion of coeff here a patch)
                    out_composition = raise_func(seq, composition)
                    if parent_basis.__class__.__name__ in (
                            'SymmetricFunctionAlgebra_homogeneous_with_category',
                            'SymmetricFunctionAlgebra_elementary_with_category',
                            'SymmetricFunctionAlgebra_power_with_category',
                            'SymmetricFunctionAlgebra_witt_with_category',
                            'SymmetricFunctionAlgebra_schur_with_category',
                            'HallLittlewood_qp_with_category'):
                        return coeff * straighten(parent_basis, out_composition)
                    else:
                        return coeff * parent_basis(out_composition)
            def call_monomial(seq, coeff, operand, power=1):
                for _ in range(power):
                    operand = raise_func(seq, operand)
                return (operand, coeff)
            # start here
            if hasattr(operand, '_get_indices_for_index_operator'):
                indices = operand._get_indices_for_index_operator()
                # TODO: see if this works for TUPLES of indices.  This has only been tested for a single index.
                new_indices = self.__call__(indices)[0][0]
                out = operand._new_object_for_index_operator(new_indices)
                return out
            elif isinstance(operand, tuple):
                # the operand is actually a tuple of operands, so perform __call__ on each piece
                return tuple(self.__call__(op) for op in operand)
            elif is_sequence(operand):
                # the operand is some kind of composition
                return [call_monomial(index, coeff, operand) for index, coeff in self]
            else:
                # the operand is a symmetric function
                if len(operand) > 1:
                    # the operand looks like s[2, 1] + s[3], for example
                    return sum(self.__call__(summand) for summand in summands(operand))
                else:
                    out_list = [call_monomial(index, coeff, operand) for index, coeff in self]
                    return sum(coeff * mon for mon, coeff in out_list)


class RaisingOperatorAlgebra(ShiftingOperatorAlgebra):
    r"""
    We follow the following convention!:

    R[(1, 0, -1)] is the raising operator that raises the first part by 1 and lowers the third part by 1.

    For a definition of raising operators, see [cat]_ Definition 2.1, but be wary that the notation is different there.  See :meth:`ij` for a way to create operators using the notation in the paper.

    If you do NOT want any restrictions on the allowed sequences, simply use 'ShiftingOperatorAlgebra' instead of 'RaisingOperatorAlgebra'.

    OPTIONAL ARGUMENTS:

    - ``base_ring`` -- (default ``QQ['t']``) the ring you will use on the raising operators.

    - ``prefix`` -- (default ``"R"``) the label for the raising operators.

    EXAMPLES::

        sage: R = RaisingOperatorAlgebra()
        sage: s = SymmetricFunctions(QQ['t']).s()
        sage: h = SymmetricFunctions(QQ['t']).h()

        sage: R[(1, -1)]
        R(1, -1)
        sage: R[(1, -1)](s[5, 4])
        s[6, 3]
        sage: R[(1, -1)](h[5, 4])
        h[6, 3]

        sage: (1 - R[(1,-1)]) * (1 - R[(0,1,-1)])
        R() - R(0, 1, -1) - R(1, -1) + R(1, 0, -1)
        sage: ((1 - R[(1,-1)]) * (1 - R[(0,1,-1)]))(s[2, 2, 1])
        (-3*t-2)*s[] + s[2, 2, 1] - s[3, 1, 1] + s[3, 2]

    """
    def __init__(self, base_ring=QQ['t'], prefix='R'):
        ShiftingOperatorAlgebra.__init__(self,
            base_ring=base_ring,
            prefix=prefix,
            basis_indices=RaisingSequenceSpace())

    def ij(self, i, j):
        r""" Shorthand element constructor that allows you to create raising operators using the familiar `R_{ij}` notation found in [cat]_ Definition 2.1, with the exception that indices here are 0-based, not 1-based.

        EXAMPLES::

            # create the raising operator which raises part 0 and lowers part 2 (indices are 0-based)
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


class PieriOperatorAlgebra(ShiftingOperatorAlgebra):
    r""" The Pieri operator `u_i`.

    EXAMPLES::

        sage: u = PieriOperatorAlgebra()

        # act on catalan function
        sage: cf = CatalanFunction([(0,2), (1,2)], [6, 6, 6])
        sage: u.i(2)(cf)
        CatalanFunction([(0,2), (1,2)], [6, 6, 5])

        # act on k schur function
        sage: base_ring = QQ['t']
        sage: Sym = SymmetricFunctions(base_ring)
        sage: t = base_ring.gen()
        sage: ks = Sym.kBoundedSubspace(4, t).kschur()
        sage: u.i(2)(ks[2, 2, 1])
        ks4[2, 2, 1] + t^2*ks4[3, 2] + t^3*ks4[4, 1]
        # TODO: verify by hand that above is really correct, or maybe a simpler example

    """
    def __init__(self, base_ring=QQ['t'], prefix='u'):
        ShiftingOperatorAlgebra.__init__(self,
            base_ring=base_ring,
            prefix=prefix,
            basis_indices=ShiftingSequenceSpace())

    def i(self, i):
        r""" Shorthand element constructor that allows you to create Pieri operators using the familiar `u_i` notation, with the exception that indices here are 0-based, not 1-based.

        EXAMPLES::

            # create the Pieri operator which lowers part 2 (indices are 0-based)
            sage: u.i(2)
            u((0, 0, -1))

        """
        if not i in NonNegativeIntegerSemiring():
            raise ValueError('i must be a natural number.  You input i = {i}.'.format(i=i))
        seq = [0] * (i + 1)
        seq[i] = -1
        seq = tuple(seq)
        return self._element_constructor_(seq)

    class Element(ShiftingOperatorAlgebra.Element):
        def __call__(self, operand):
            if is_k_schur(operand):
                # convert to catalans
                kschur = operand.parent()
                base_ring = kschur.base_ring()
                cat_coeff_pairs = [
                    (CatalanFunctions().init_from_k_schur(kschur(index), base_ring=base_ring), coeff)
                    for index, coeff in operand]
                # act
                new_cat_coeff_pairs = [(self.__call__(cat), coeff) for cat, coeff in cat_coeff_pairs]
                # convert back to kschur
                out_coeff_pairs = [(kschur(cat.eval()), coeff) for cat, coeff in cat_coeff_pairs]
                return sum(coeff * func for func, coeff in out_coeff_pairs)
            else:
                return ShiftingOperatorAlgebra.Element.__call__(self, operand)


class HallLittlewoodVertexOperator:
    r"""
    Garsia's version of Jing's Hall-Littlewood vertex operators.  These are defined in equations 4.2 and 4.3 of [cat]_ and appear visually as a bold capital H.

    INPUTS:

    - ``base_ring`` -- (default ``QQ['t']``) the base ring to build the SymmetricFunctions upon.

    EXAMPLES::

        sage: H = HallLittlewoodVertexOperator
        sage: one = SymmetricFunctions(QQ['t']).hall_littlewood().Qp().one()
        sage: H([4, 1, 3])(one) == H(4)(H(1)(H(3)(one)))
        True

    """
    def __init__(self, composition, base_ring=QQ['t']):
        if composition in NonNegativeIntegerSemiring():
            self.composition = [composition]
        elif is_sequence(composition):
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
    r"""
    Given gamma, returns the compositional Hall-Littlewood polynomial `H_{\gamma}(\mathbf{x}; t)` in the Q' basis, as defined in [cat]_ section 4.4.

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

def raising_roots_operator(roots, t=1, base_ring=QQ['t']):
    r""" Given a list of roots `roots = \Phi` (often a root ideal), and optionally a variable `t`, return the operator

    ..  MATH::

        \prod_{(i,j) \in \Phi} (1 - tR_{ij}).

    If you input an integer for roots (e.g. ``roots = 3``), it will use the biggest possible root ideal in the `n` x `n` grid (the '`n`-th staircase root ideal').
    """
    R = RaisingOperatorAlgebra(base_ring=base_ring)
    def prod(iterable):
        return reduce(operator.mul, iterable, R.one())
    if roots in NonNegativeIntegerSemiring():
        roots = staircase_root_ideal(roots)
    op = prod([1 - t*R.ij(i, j) for (i, j) in roots])
    return op

def qt_raising_roots_operator(roots, t=None, q=None, base_ring=QQ['t', 'q']):
    r""" The q-t analogue of :meth:`raising_roots_operator`, defined by

    ..  MATH::
        \prod_{ij \in \Phi} (1 - tR_{ij}) \prod_{ij \in \Phi} (1 - qR_{ij}).

    """
    if roots in NonNegativeIntegerSemiring():
        roots = staircase_root_ideal(roots)
    if t is None:
        t = base_ring.gens()[0]
    if q is None:
        q = base_ring.gens()[1]
    op1 = raising_roots_operator(roots, t=q, base_ring=base_ring)
    op2 = raising_roots_operator(roots, t=t, base_ring=base_ring)
    return lambda x: (op2 * op1)(x)


class CatalanFunction:
    r""" A catalan function `H(\Psi; \gamma)` as discussed in [cat]_ section 4.

    The parameters for initializing a catalan function are exactly the same as in :meth:`CatalanFunctions.init_from_indexed_root_ideal`.

    EXAMPLES::

        sage: CatalanFunction([(0,2), (1,2)], [6, 6, 5])
        H([(0, 2), (1, 2)], [6, 6, 5])

    There are in fact many ways to initialize a catalan function, and the methods for doing so are found in :class:`CatalanFunctions`.

    """
    BASE_RING_DEFAULT = QQ['t']
    PREFIX_DEFAULT = 'H'

    def __init__(self, roots, index, base_ring=None, prefix=None):
        assert is_roots(roots)
        self.roots = roots
        assert is_sequence(index)
        self.index = index
        self.base_ring = base_ring if base_ring is not None else self.BASE_RING_DEFAULT
        self.prefix = prefix if prefix is not None else self.PREFIX_DEFAULT

    def __eq__(self, other):
        # TODO: account for the fact that DIFFERENT root/index pairs could actually give the SAME catalan function!!
        return set(self.roots) == set(other.roots) and self.index == other.index and self.base_ring == other.base_ring

    def __repr__(self):
        return '{}({}, {})'.format(self.prefix, self.roots, self.index)

    def _get_indices_for_index_operator(self):
        return self.index

    def _new_object_for_index_operator(self, new_index):
        new_obj = self.__class__(self.roots, new_index, base_ring=self.base_ring)
        return new_obj

    def eval(self):
        # setup
        hl = SymmetricFunctions(self.base_ring).hall_littlewood().Qp()
        t = self.base_ring.gen()
        # formula
        n = len(self.index)
        roots_complement = RI.complement(self.roots, n)
        op = raising_roots_operator(roots_complement, t=t, base_ring=self.base_ring)
        hl_poly = hl(self.index)
        cat_func = op(hl_poly)
        return cat_func


class CatalanFunctions:
    r""" The family of Catalan Functions.

    Use this class as a factory to initialize a :class:`CatalanFunction` object with any valid identifying data.  See the ``init_from...`` methods below for all possible ways to create a catalan function.
    """

    def init_from_indexed_root_ideal(self, roots, index, base_ring=None, prefix=None):
        r"""
        INPUTS:

        - ``roots`` -- iterable of roots `\Phi` (typically a root ideal)

        - ``index`` -- composition `\gamma` that indexes the root ideal and appears in the Hall-Littlewood Q' function `H_\gamma`

        OPTIONAL INPUTS:

        - ``base_ring`` -- (default ``QQ['t']``) the ring over which to build the `h_\gamma(x; \alpha)`'s

        OUTPUTS:

        The catalan function

        ..  MATH::

            \prod_{ij \in \Phi} (1 - R_{ij}) H_\gamma

        EXAMPLES::

            sage: CatalanFunctions().init_from_indexed_root_ideal([(0,2), (1,2)], [6, 6, 5])
            H([(0, 2), (1, 2)], [6, 6, 5])

        """
        return CatalanFunction(roots, index, base_ring, prefix)

    def init_from_skew_partition(self, sp, base_ring=None, prefix=None):
        r""" Given a SkewPartition ``sp``, return the catalan function `H(\Phi^+(sp); rs(sp))`.
        """
        ri = skew_partition_to_root_ideal(sp, type='max')
        rs = sp.row_lengths()
        return self.init_from_indexed_root_ideal(ri, rs, base_ring, prefix)

    def init_from_row_and_column_lengths(self, row_lengths, column_lengths, base_ring=None, prefix=None):
        r""" Determine the skew partition `D` with row-shape ``row_lengths`` and column-shape ``column_lengths``, and return the catalan function `H(\Phi^+(D); \text{row_lengths})`.
        """
        sp = SkewPartitions().init_from_row_and_column_length(row_lengths, column_lengths)
        return self.init_from_skew_partition(sp, base_ring, prefix)

    def init_from_k_shape(self, p, k, base_ring=None, prefix=None):
        r""" Given `k` and a `k`-shape `p`, return the catalan function `H(\Phi^+(p); rs(p))`.
        """
        assert is_k_shape(p, k)
        sp = SkewPartition([p, []])
        return self.init_from_skew_partition(self, sp, base_ring, prefix)

    def init_from_k_schur(self, func, base_ring=None, prefix=None):
        r""" Given a k-schur function ``func`` `= s^k_\lambda(x;t)`, initialize the catalan function `H(\Delta^k(\lambda); \lambda)`.

        Mathematically, these two functions are equal.  The usefulness of this method is that you input a ``sage.combinat.sf.new_kschur.kSchur_with_category`` object and you obtain a ``CatalanFunction`` object.

        EXAMPLES::

            sage: base_ring = QQ['t']
            sage: Sym = SymmetricFunctions(base_ring)
            sage: t = base_ring.gen()
            sage: ks = Sym.kBoundedSubspace(4, t).kschur()
            sage: func = ks[2, 1, 1]
            sage: CatalanFunction(func, base_ring=base_ring)
            H([], [2, 1, 1])

            # TODO: make sure [] above is correct.  go by hand.
        """
        # check inputs
        assert is_k_schur(func)
        assert len(func.support()) == 1
        # gather roots
        index = func.support()[0]
        k = func.parent().k
        roots = partition_to_k_schur_root_ideal(index, k)
        # return
        return self.init_from_indexed_root_ideal(roots, index, base_ring, prefix)


##############
def k_plus_one_core_to_k_schur_function(p, k, base_ring=QQ['t']):
    # TODO: compare the performance of this function to existing k-schur function.
    p = Core(p, k+1)
    return k_shape_to_catalan_function(p, k, base_ring)


class InfiniteDimensionalFreeAlgebra(CombinatorialFreeModule):
    r"""
    By default, the algebra generated by ``x[0], x[1], x[2], ...`` over the integers.

    To change the index set of the generators, use ``index_set=`` (default ``NN``).  To overhaul the set of generators entirely (not recommended), use ``basis_indices=``.

    To change the ring that the algebra works over, use ``base_ring=`` (default ``ZZ``).

    To change the prefix of the generators, use ``prefix=`` (default ``'x'``).
    """
    def __init__(self,
            base_ring=IntegerRing(),
            prefix='x',
            basis_indices=None,
            index_set=NonNegativeIntegerSemiring()):
        self._base_ring = base_ring
        self._basis_monoid = FreeMonoid(index_set=index_set, commutative=True, prefix=prefix) if basis_indices is None else basis_indices
        # category
        category = Algebras(self._base_ring.category()).WithBasis().Commutative()
        category = category.or_subcategory(category)
        # init
        CombinatorialFreeModule.__init__(
            self,
            self._base_ring,
            self._basis_monoid,
            category=category,
            prefix='',
            bracket=False)
        # TODO: make a SEPARATE class called InfiniteDimensionalFreeRing or similar
        # self._init_category_(CommutativeRings()) # i think .Commutative() above is a better solution

    def is_prime_field(self):
        return False

    def _element_constructor_(self, monoid_el):
        assert monoid_el in self._basis_monoid
        return self.basis()[monoid_el]

    def __getitem__(self, user_input):
        # USER front entrance to creating elements "x[4]"
        assert user_input in IntegerRing()
        monoid_el = self._basis_monoid.gen(user_input)
        return self.basis()[monoid_el]

    @cached_method
    def one_basis(self):
        # identity index
        return self._basis_monoid.one()

    def product_on_basis(self, monoid_el1, monoid_el2):
        monoid_el_product = monoid_el1 * monoid_el2
        return self._element_constructor_(monoid_el_product)

    def _repr_(self):
        return "{class_name} with generators indexed by integers, over {base_ring}".format(class_name=self.__class__.__name__, base_ring=self._base_ring)


DoubleRing = InfiniteDimensionalFreeAlgebra(prefix='a', index_set=IntegerRing())
r"""   ``DoubleRing`` is the ring `\Lambda(a)` found in [Fun]_ section 3. """


def dual_k_theoretic_homogeneous(k, r, base_ring=QQ):
    r""" The dual K-theoretic h, often denoted Kh, is defined for any integer `k` by the formula `h_k(x, r) = \sum_{i=0}^{k} \binom{r + i - 1}{i} h_{k - i}(x)` in [LN]_ p.88 top-right.

    If `k` and `r` are compositions, then it is recursively defined as `h_k(x, r) = \prod_j h_{k_j}(x, r_j)`.

    EXAMPLES::

        sage: dual_k_theoretic_homogeneous(0, 0)
        1

        sage: dual_k_theoretic_homogeneous(1, 2, base_ring=QQ['t'])
        h[1] + 2

        sage: dual_k_theoretic_homogeneous([2, 1], [1, 1])
        h[1]**2 + h[1]*h[2] + 2*h[1] + h[2] + 1

    """
    if is_sequence(k):
        # pad with 0's
        max_len = max(len(k), len(r))
        k = list(k) + [0] * (max_len - len(k))
        r = list(r) + [0] * (max_len - len(r))
        # multiply
        h_list = [dual_k_theoretic_homogeneous(k_el, r_el, base_ring) for k_el, r_el in zip(k, r)]
        return reduce(operator.mul, h_list)
    else:
        assert k >= 0
        h = SymmetricFunctions(base_ring).h()
        return sum(binomial(r + i - 1, i) * h[k - i] for i in range(k + 1))

def dual_k_catalan_function(roots, index, index2, base_ring=QQ):
    r"""
    INPUTS:

    - ``roots`` -- iterable of roots `\Phi` (typically a root ideal)

    - ``index`` -- composition `\gamma` that indexes `h_\gamma(x; \alpha)`

    - ``index2`` -- composition `\alpha` used in `h_\gamma(x; \alpha)`

    OPTIONAL INPUTS:

    - ``base_ring`` -- (default ``QQ``) the ring over which to build the `h_\gamma(x; \alpha)`'s

    OUTPUTS:

    The 'dual k' catalan function

    ..  MATH::

        \prod_{ij \in \Delta^+ \smallsetminus \Phi} (1 - R_{ij}) h_\gamma(x; \alpha).

    """
    # setup
    Kh = dual_k_theoretic_homogeneous(index, index2, base_ring=base_ring)
    # formula
    n = len(index)
    roots_complement = RI.complement(roots, n)
    op = raising_roots_operator(roots_complement, t=1, base_ring=base_ring)
    cat_func = op(Kh)
    return cat_func

def dual_grothendieck_function(composition, base_ring=QQ):
    r""" Given a composition `composition = \lambda`, return the dual Grothendieck function defined by `g_\lambda(x) = \text{det}(h_{\lambda_i + j - i}(x, i - 1))` in [LN]_ p.88 equation (4).

    Equivalently, the dual Grothendieck function is `g_\lambda(x) = \prod_{ij \in \Delta^+} (1 - R_{ij}) h_\lambda(x; (0, 1, \ldots, n-1))`.

    EXAMPLES::

        sage: h = SymmetricFunctions(QQ).h()
        sage: dual_grothendieck_function([2, 1])
        h[1]*h[2] + h[2] - h[3]

    """
    roots = [] # because dual_k_catalan_function will take the complement
    n = len(composition)
    reversed_staircase_ptn = list(reversed(staircase_shape(n)))
    return dual_k_catalan_function(roots, composition, reversed_staircase_ptn, base_ring=base_ring)

def double_homogeneous_building_block(p, n):
    r""" The double complete homogeneous symmetric polynomial building block `h_p(x \,||\, a)` defined as

    ..  MATH::
        h_p(x_1, \ldots, x_n \,||\, a) \,= \sum_{n \geq i_1 \geq \ldots \geq i_p \geq 1} (x_{i_1} - a_{i_1})(x_{i_2} - a_{i_2 - 1}) \cdots (x_{i_p} - a_{i_p - p + 1})

    in [Fun]_ section 3 between equation (6) and (7).  Note that our indices are 0-based.
    """
    def prod(iterable):
        return reduce(operator.mul, iterable, a.one())
    a = DoubleRing
    sym = SymmetricFunctions(DoubleRing)
    s = sym.s()
    one = s.one()
    one_poly = one.expand(n)
    x = one_poly.parent().gens()
    ptns = Partitions(max_part=n, length=p)
    total_sum = 0
    for ptn in ptns:
        summand = prod([(x[ptn[b]] - a[ptn[b] - b]) for b in range(p)])
        total_sum += summand
    return total_sum

def shift(element):
    r""" The function `\tau` which acts on any element of `\Lambda(a)` (``DoubleRing``) by sending each element `a_i` to `a_{i+1}` for all `i`.  It can be found in [Fun]_ p.8 between equations (6) and (7).
    """
    # idea, use .hom ``f = ZZ.hom(GF(3))``
    # sage: R.<x> = ZZ[]
    # sage: f = R.hom([x])

    # idea, define something on x[i] in general and take the induces hom

    # idea, manual
    a = DoubleRing
    new_monomials = []
    for monomial, coeff in element.monomialomial_coefficients():
        new_factors = []
        for factor in monomial:
            new_factor = a(factor.index() + 1)
            new_factors.append(new_factor)
        new_monomial = coeff * prod(new_factors)
        new_monomials.append(new_monomial)
    new_element = a.sum(new_monomials)
    return new_element

def double_homogeneous_building_block_shifted(r, s, n):
    r""" Given `r` and `s`, returns `h_{r, s} = \tau^s h_r(x \,||\, a)`, as defined in [Fun]_ before eq (8).
    """
    if n > 0:
        out = double_homogeneous_building_block(r, n)
        for _ in range(s):
            out = shift(out)
        return out
    else:
        raise NotImplemented

class DoubleHomogeneous:
    r"""
    mu1 -- composition
    mu2 -- composition
    n -- number of `x` variables

    EXAMPLE::
        # create the double homogeneous `h^{(4)}_{\mu, \beta}`
        sage: DoubleHomogeneous(mu, beta, 4)

        # create the double homogeneous shifted building block `h_{r, s}` in 4 variables
        sage: r = 5
        sage: s = 2
        sage: DoubleHomogeneous([r], [s], 4)

        # create the double homogeneous symmetric building block `h_p(x \,||\, a)` in 4 variables
        sage: p = 3
        sage: DoubleHomogeneous([p], [0], 4)
    """
    def __init__(self, index1, index2, prefix='h'):
        self.index1 = index1
        self.index2 = index2
        self.prefix = prefix

    def __repr__(self):
        r""" For example, h(n)[mu, beta], which represents `h^{(n)}_{\mu, \beta}`. """
        return '{}({})[{}, {}]'.format(self.prefix, self.n, self.index1, self.index2)

    def _get_indices_for_index_operator(self):
        pair = (self.index1, self.index2)
        return pair

    def _new_object_for_index_operator(self, indices):
        (index1, index2) = indices
        new_obj = self.__class__(index1, index2, self.n)
        return new_obj

    def eval(self, n):
        r""" Given the number of variables `n`, output self expanded in terms of the shifted double homogeneous building blocks `h_{r, s}`. """
        (mu1, mu2) = (self.index1, self.index2)
        # pad with 0's
        max_len = max(len(mu1), len(mu2))
        mu1 = list(mu1) + [0] * (max_len - len(mu1))
        mu2 = list(mu2) + [0] * (max_len - len(mu2))
        # compute
        h_list = [double_homogeneous_shifted(mu1[i], mu2[i], n) for i in range(max_len)]
        hp = prod(h_list)
        return hp

def double_schur(index, n):
    r""" Given a composition `index = \lambda` and the number of variables `n`, return the double Schur function defined by

    ..  MATH::
        s_\lambda(x_1, \ldots, x_n \,||\, a) = \text{det}\left(h^{(n)}_{\lambda_i + i - j, j - 1}\right)

    or equivalently, defined by

    ..  MATH::
        s_\lambda(x_1, \ldots, x_n \,||\, a) = \prod_{ij \in \Delta^+(l)} (1 - R_{ij}) h^{(n)}_{\lambda, (0, \ldots, l-1)}

    where `l` is the length of `\lambda`, in [Fun]_ p.9 equation (9).

    """
    l = len(index)
    op = raising_roots_operator(l, t=1, base_ring=ZZ)
    rho = list(reversed(staircase_shape(l)))
    h_index = DoubleHomogeneous(index, rho, n)
    return op(h_index)

def double_catalan_function(roots, index, n):
    # setup
    l = len(index)
    rho = list(reversed(staircase_shape(l)))
    h_index = DoubleHomogeneous(index, rho, n)
    # formula
    roots_complement = RI.complement(roots, l)
    # TODO: see what to pass in for base_ring in below line.
    op = raising_roots_operator(roots_complement, t=1)
    cat_func = op(h_index)
    return cat_func
