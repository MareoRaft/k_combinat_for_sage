# -*- coding: utf-8 -*-
r"""
This module contains all functionalities that are not already organized into the other files.  New functionalities written to the library often appear here, and eventually get organized into separate files.

AUTHORS:

- Matthew Lancellotti (2018): Initial version

REFERENCES:

.. [Fun] `Raising operators and the Littlewood-Richardson polynomials <https://arxiv.org/pdf/1203.4729.pdf>`_.  Fun, Alex.
.. [LN] `Finite sum Cauchy identity for dual Grothendieck polynomials <https://projecteuclid.org/download/pdf_1/euclid.pja/1407415930>`_.
"""

#*****************************************************************************
#  Copyright (C) 2018 Matthew Lancellotti <mvlancellotti@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.all import *

from core import *
import core
from partition import *
import partition
from partition import _is_sequence
from skew_partition import *
import skew_partition
from k_shape import *
import k_shape
from root_ideal import *
import root_ideal
from strong_marked_tableau import *
import strong_marked_tableau
# ^*^ sphinx insert ^*^


# HELPERS
def _is_k_schur(obj):
    r""" Helper function.

    Checks if ``obj`` is coming from the 'kSchur_with_category' class.

    EXAMPLES::

    An example of a `k`-schur function::

        sage: Sym = SymmetricFunctions(QQ)
        sage: KB = Sym.kBoundedSubspace(3,1)
        sage: ks = KB.kschur()
        sage: _is_k_schur(ks[2, 1] + ks[1, 1])
        True

    The following is a schur function, *not* a `k`-schur function:

        sage: Sym = SymmetricFunctions(QQ)
        sage: s = Sym.schur()
        sage: _is_k_schur(s[2, 1] + s[1, 1])
        False
    """
    try:
        classname = obj.parent().__class__.__name__
        return classname == 'kSchur_with_category'
    except:
        return False


# MAIN
def get_k_rectangles(k):
    r""" Return the list of ``k``-rectangles.

    A __``k``-rectangle__ is a partition whose Ferrer's diagram is a rectangle whose largest hook-length is `k`.

    EXAMPLES::

        sage: get_k_rectangles(0)
        []
        sage: get_k_rectangles(1)
        [[1]]
        sage: get_k_rectangles(2)
        [[2], [1, 1]]
        sage: get_k_rectangles(3)
        [[3], [2, 2], [1, 1, 1]]
    """
    return [Partition([a] * b) for (a, b) in k_rectangle_dimension_list(k)]


def get_k_irreducible_partition_lists(k):
    r""" Return the list of ``k``-irreducible partitions.

    The `k`-irreducible partitions are output at lists, not Partition objects.

    There are `k!` such partitions, and computation time starts to get slow around `k = 10`.

    EXAMPLES::

        sage: get_k_irreducible_partition_lists(0)
        [[]]
        sage: get_k_irreducible_partition_lists(1)
        [[]]
        sage: get_k_irreducible_partition_lists(2)
        [[], [1]]
        sage: get_k_irreducible_partition_lists(3)
        [[], [1], [1, 1], [2], [2, 1], [2, 1, 1]]

    ..  SEEALSO::

        :meth:`get_k_irreducible_partitions`
    """
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
    r""" Return the list of ``k``-irreducible partitions.

    The output is a list of :class:`Partition` objects.

    There are `k!` such partitions, and computation time starts to get slow around `k = 10`.

    EXAMPLES::

        sage: get_k_irreducible_partitions(0)
        [[]]
        sage: get_k_irreducible_partitions(1)
        [[]]
        sage: get_k_irreducible_partitions(2)
        [[], [1]]
        sage: get_k_irreducible_partitions(3)
        [[], [1], [1, 1], [2], [2, 1], [2, 1, 1]]

    ..  SEEALSO::

        :meth:`get_k_irreducible_partition_lists`
    """
    return [Partition(e) for e in get_k_irreducible_partition_lists(k)]


def size_to_num_linked_partition_self_pairs(size):
    r""" Given a natural number ``size``, count how many partitions `l` of size ``size`` have the property that `(l, l)` has a corresponding skew-linked-diagram.

    Note: A 'skew-linked-diagram' is a :class:`SkewPartition` that is linked.

    EXAMPLES::

        sage: size_to_num_linked_partition_self_pairs(0)
        1
        sage: size_to_num_linked_partition_self_pairs(1)
        1
        sage: size_to_num_linked_partition_self_pairs(2)
        1
        sage: size_to_num_linked_partition_self_pairs(3)
        2

    ..  SEEALSO::

        :meth:`SkewPartition.is_linked`
    """
    # DO NOT ADD TO SAGE
    ps = Partitions(size)
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
    # DO NOT ADD TO SAGE
    n = 0
    while n < num_terms:
        print('n={}\t{}=f(n)'.format(n, func(n)))


def size_to_k_shapes(n, k):
    # DO NOT ADD TO SAGE
    r""" Return all partitions of size ``n`` that are ``k``-shapes. """
    return [ptn for ptn in Partitions(n) if is_k_shape(ptn, k)]


def size_to_num_k_shapes(n, k):
    # DO NOT ADD TO SAGE
    return len(size_to_k_shapes(n, k))


def straighten(basis, gamma):
    r""" Perform Schur function straightening by the Schur straightening rule.

    See [cat]_, Prop. 4.1.  Also known as the slinky rule.

    .. MATH::

        s_{\gamma}(\mathbf{x}) = \begin{cases}
            \text{sgn}(\gamma+\rho) s_{\text{sort}(\gamma+\rho) -\rho}(\mathbf{x}) & \text{if } \gamma + \rho \text{ has distinct nonnegative parts,} \\
            0 & \text{otherwise,}
        \end{cases}

    where `\rho=(\ell-1,\ell-2,\dots,0)`, `\text{sort}(\beta)` denotes the weakly decreasing sequence obtained by sorting `\beta`, and `\text{sgn}(\beta)` denotes the sign of the (shortest possible) sorting permutation.

    EXAMPLES:

    We know s[2, 1, 3] := -s[2, 2, 2]::

        sage: s = SymmetricFunctions(QQ).s()
        sage: straighten(s, [2, 1, 3])
        -s[2, 2, 2]
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
    gamma = Composition(gamma)
    if basis.__class__.__name__ in ('SymmetricFunctionAlgebra_monomial_with_category', 'SymmetricFunctionAlgebra_dual_with_category'):
        raise NotImplemented(
            'Straightening does not exist (that i know of) for the monomial basis or the forgotten/dual basis.')
    elif basis.__class__.__name__ == 'HallLittlewood_qp_with_category':
        return compositional_hall_littlewood_Qp(gamma, base_ring=basis.base_ring(), t=basis.t)
    elif basis.__class__.__name__ in ('SymmetricFunctionAlgebra_homogeneous_with_category', 'SymmetricFunctionAlgebra_elementary_with_category', 'SymmetricFunctionAlgebra_power_with_category', 'SymmetricFunctionAlgebra_witt_with_category'):
        new_gamma = list(reversed(sorted(gamma)))
        if has_nonnegative_parts(new_gamma):
            return basis(Partition(new_gamma))
        else:
            return 0
    elif basis.__class__.__name__ == 'SymmetricFunctionAlgebra_schur_with_category':
        rho = list(range(len(gamma) - 1, -1, -1))
        combined = [g + r for g, r in zip(gamma, rho)]
        if has_distinct_parts(combined) and has_nonnegative_parts(combined):
            sign = (-1)**number_of_noninversions(combined)
            sort_combined = reversed(sorted(combined))
            new_gamma = [sc - r for sc, r in zip(sort_combined, rho)]
            return sign * basis(Partition(new_gamma))
        else:
            return 0
    else:
        raise ValueError(
            "The input parameter 'basis' should be a symmetric function basis.  For example, 's = SymmetricFunctions(QQ).s(); straighten(s, [2, 1, 3])', or one could use 'h' instead of 's'.")


class ShiftingSequenceSpace():
    r""" A helper for :class:`ShiftingOperatorAlgebra.`

    Helps :class:`ShiftingOperatorAlgebra` know which indices are valid and which indices are not for the basis.

    EXAMPLES::

        sage: S = ShiftingSequenceSpace()
        sage: (1, -1) in S
        True
        sage: (1, -1, 0, 9) in S
        True
        sage: [1, -1] in S
        False
        sage: (0.5, 1) in S
        False
    """
    def __init__(self, base=IntegerRing()):
        self.base = base
        # category = InfiniteEnumeratedSets()
        # Parent.__init__(self, category=category)

    def __contains__(self, seq):
        r""" Returns ``True`` if and only if ``seq`` is a valid shifting sequence.

        EXAMPLES::

            sage: S = ShiftingSequenceSpace()
            sage: (1, -1) in S
            True
            sage: (1, -1, 0, 9) in S
            True
            sage: [1, -1] in S
            False
            sage: (0.5, 1) in S
            False
        """
        if not isinstance(seq, tuple):
            return False
        return not any(i not in self.base for i in seq)

    CHECK_ERROR_MESSAGE = 'Expected valid index (a tuple of {base}), but instead received {seq}.'

    def check(self, seq):
        r""" Verify that ``seq`` is a valid shifting sequence.

        If it is not, raise an error.

        EXAMPLES::

            sage: S = ShiftingSequenceSpace()
            sage: S.check((1, -1))
            sage: S.check((1, -1, 0, 9))
            sage: S.check([1, -1])
            Traceback (most recent call last):
            ...
            ValueError: Expected valid index (a tuple of Integer Ring), but instead received [1, -1].
            sage: S.check((0.5, 1))
            Traceback (most recent call last):
            ...
            ValueError: Expected valid index (a tuple of Integer Ring), but instead received [1, -1].
        """
        if not self.__contains__(seq):
            raise ValueError(self.CHECK_ERROR_MESSAGE.format(
                base=self.base, seq=seq))


class RaisingSequenceSpace(ShiftingSequenceSpace):
    r""" A helper for :class:`RaisingOperatorAlgebra`.

    Helps :class:`RaisingOperatorAlgebra` know which indices are valid and which indices are not for the basis.

    EXAMPLES::

        sage: RS = RaisingSequenceSpace()
        sage: (1, -1) in RS
        True
        sage: (1, 0, -1) in RS
        True
        sage: (1, -1, 0, 9) in RS
        False
        sage: [1, -1] in RS
        False
    """

    CHECK_ERROR_MESSAGE = 'Expected valid index (a tuple of {base} elements, where every partial sum is nonnegative and every total sum is 0), but instead received {seq}.'

    def __contains__(self, seq):
        r""" Returns ``True`` if and only if ``seq`` is a valid raising sequence.

        EXAMPLES::

            sage: RS = RaisingSequenceSpace()
            sage: (1, -1) in RS
            True
            sage: (1, 0, -1) in RS
            True
            sage: (1, -1, 0, 9) in RS
            False
            sage: [1, -1] in RS
            False
        """
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
    r""" An algebra of shifting operators.

    We follow the following convention:

    S[(1, 0, -1, 2)] is the shifting operator that raises the first part by 1, lowers the third part by 1, and raises the fourth part by 2.

    OPTIONAL ARGUMENTS:

    - ``base_ring`` -- (default ``QQ['t']``) the ring you will use on the raising operators.

    - ``prefix`` -- (default ``"R"``) the label for the raising operators.

    EXAMPLES::

        sage: S = ShiftingOperatorAlgebra()
        sage: s = SymmetricFunctions(QQ['t']).s()
        sage: h = SymmetricFunctions(QQ['t']).h()

        sage: S[(1, -1, 2)]
        S(1, -1, 2)
        sage: S[(1, -1, 2)](s[5, 4])
        s[6, 3, 2]
        sage: S[(1, -1, 2)](h[5, 4])
        h[6, 3, 2]

        sage: (1 - S[(1,-1)]) * (1 - S[(4,)])
        S() - S(1, -1) - S(4,) + S(5, -1)
        sage: ((1 - S[(1,-1)]) * (1 - S[(4,)]))(s[2, 2, 1])
        s[2, 2, 1] - s[3, 1, 1] - s[6, 2, 1] + s[7, 1, 1]

    ..  SEEALSO::
        :class:`RaisingOperatorAlgebra`, :class:`PieriOperatorAlgebra`
    """
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
        r""" Return the shifting operator whose index is ``seq``.

        This method is only for basis indices.

        EXAMPLES::

            sage: S = ShiftingOperatorAlgebra()
            sage: S[1, 1, -9]
            S(1, 1, -9)
        """
        # seq should be a basis index
        self._basis_indices.check(seq)
        return self.basis()[seq]

    def _element_constructor_(self, seq):
        r""" Return the shifting operator whose index is ``seq``.

        This method is only for basis indices.

        EXAMPLES::

            sage: S = ShiftingOperatorAlgebra()
            sage: S._element_constructor_([1, 1, -9])
            S(1, 1, -9)
            sage: S[1, 1, -9]
            S(1, 1, -9)
        """
        return self.__getitem__(seq)

    @cached_method
    def one_basis(self):
        r""" Return the index of the identity element.

        EXAMPLES::

            sage: S = ShiftingOperatorAlgebra()
            sage: S.one_basis()
            ()
        """
        return tuple()

    def _repr_(self):
        r""" Return a string describing ``self`` to humans.

        EXAMPLES::

            sage: S = ShiftingOperatorAlgebra()
            sage: S
            Shifting Operator Algebra over Univariate Polynomial Ring in t over Rational Field
        """
        return "Shifting Operator Algebra over {base_ring}".format(base_ring=self._base_ring)

    def product_on_basis(self, index1, index2):
        r""" Given indices ``index1`` and ``index2``, return the product of the basis elements indexed by those indices.

        EXAMPLES::

            sage: S = ShiftingOperatorAlgebra()
            sage: S.product_on_basis([1, 1], [0, 1, 2])
            S(1, 2, 2)
        """
        # pad with 0's
        max_len = max(len(index1), len(index2))
        index1 = tuple(index1) + (0,) * (max_len - len(index1))
        index2 = tuple(index2) + (0,) * (max_len - len(index2))
        # add the vectors
        index_product = tuple(i1 + i2 for i1, i2 in zip(index1, index2))
        return self.__getitem__(index_product)

    class Element(CombinatorialFreeModule.Element):
        r""" An element of a :class`ShiftingOperatorAlgebra`. """

        def indices(self):
            r""" Return the support of ``self``.

            EXAMPLES::

                sage: S = ShiftingOperatorAlgebra()
                sage: (S[2, 1] + S[1, 1]).indices()
                [(1, 1), (2, 1)]
            """
            return self.support()

        def index(self):
            r""" Return the index of ``self``.

            This method is only for basis elements.

            EXAMPLES::

                sage: S = ShiftingOperatorAlgebra()
                sage: S[2, 1].index()
                (2, 1)
                sage: (S[2, 1] + S[1, 1]).index()
                Traceback (most recent call last):
                ...
                ValueError: This is only defined for basis elements.  For other elements, use indices() instead.
         """
            if len(self) != 1:
                raise ValueError(
                    "This is only defined for basis elements.  For other elements, use indices() instead.")
            return self.indices()[0]

        def _call_basis_on_index(self, seq, index):
            """ For internal use only!

            Return the action of the basis element indexed by ``seq`` upon the composition ``index``.

            INPUTS:

            - ``seq`` -- The sequence of the basis element that acts.

            - ``index`` -- A sequence (typically a composition or a partition) that we act upon.

            EXAMPLES::

                sage: S = ShiftingOperatorAlgebra()
                sage: S[2, 1]._call_basis_on_index([1, 1], [1, 2, 3, 4, 5])
                [2, 3, 3, 4, 5]

            ..  SEEALSO::
                :meth:`__call__`
            """
            assert _is_sequence(index)
            # pad sequence and index with 0's
            index = list(index) + [0] * (len(seq) - len(index))
            seq = tuple(seq) + (0,) * (len(index) - len(seq))
            # raise and drop
            return [v + s for v, s in zip(index, seq)]

        def __call__(self, operand):
            r""" Return the action of this shifting operator element on the index ``operand``.

            EXAMPLES::

                sage: S = ShiftingOperatorAlgebra()

                sage: S[2, 1]([1, 2, 3, 4])
                [([3, 3, 3, 4], 1)]

                sage: S[2, 1](([1, 2, 3, 4], [1, 1, 1]))
                ([([3, 3, 3, 4], 1)], [([3, 2, 1], 1)])

                sage: Sym = SymmetricFunctions(QQ['t'])
                sage: s = Sym.schur()
                sage: S[2, 1](s[4, 3, 2, 1])
                s[6, 4, 2, 1]
            """
            def raise_func(seq, operand):
                # seq is the index
                if _is_sequence(operand):
                    return self._call_basis_on_index(seq, operand)
                else:
                    # it's some symmetric function basis element
                    parent_basis = operand.parent()
                    # but does the base ring of the parent basis agree with the base ring of the operator algebra??
                    # process the vectors
                    dic = operand.monomial_coefficients()
                    assert len(dic) == 1
                    # occasionally a coefficient can show up (not cool, so consider the inclusion of coeff here a patch)
                    (composition, coeff) = dic.items()[0]
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
            elif _is_sequence(operand):
                # the operand is some kind of composition
                return [call_monomial(index, coeff, operand) for index, coeff in self]
            else:
                # the operand is a symmetric function
                if len(operand) > 1:
                    # the operand looks like s[2, 1] + s[3], for example
                    return sum(self.__call__(summand) for summand in operand.terms())
                else:
                    out_list = [call_monomial(index, coeff, operand)
                                for index, coeff in self]
                    return sum(coeff * mon for mon, coeff in out_list)


class RaisingOperatorAlgebra(ShiftingOperatorAlgebra):
    r""" An algebra of raising operators.

    This class subclasses :class:`ShiftingOperatorAlgebra` and inherits the large majority of its functionality from there.

    We follow the following convention!:

    R[(1, 0, -1)] is the raising operator that raises the first part by 1 and lowers the third part by 1.

    For a definition of raising operators, see [cat]_ Definition 2.1, but be wary that the notation is different there.  See :meth:`ij` for a way to create operators using the notation in the paper.

    If you do NOT want any restrictions on the allowed sequences, use :class:`ShiftingOperatorAlgebra` instead of :class:`RaisingOperatorAlgebra`.

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

    ..  SEEALSO::
        :class:`ShiftingOperatorAlgebra`
    """

    def __init__(self, base_ring=QQ['t'], prefix='R'):
        ShiftingOperatorAlgebra.__init__(self,
                                         base_ring=base_ring,
                                         prefix=prefix,
                                         basis_indices=RaisingSequenceSpace())

    def ij(self, i, j):
        r""" Return the raising operator `R_{ij}` as notated in [cat]_ Definition 2.1.

        Shorthand element constructor that allows you to create raising operators using the familiar `R_{ij}` notation found in [cat]_ Definition 2.1, with the exception that indices here are 0-based, not 1-based.

        EXAMPLES:

        Create the raising operator which raises part 0 and lowers part 2 (indices are 0-based)::

            sage: R.ij(0, 2)
            R((1, 0, -1))

        ..  SEEALSO::
            :meth:`ShiftingOperatorAlgebra._element_constructor_`, :meth:`ShiftingOperatorAlgebra.__getitem__`
        """
        if not i in NonNegativeIntegerSemiring():
            raise ValueError(
                'i must be a natural number.  You input i = {i}.'.format(i=i))
        if not j in NonNegativeIntegerSemiring():
            raise ValueError(
                'j must be a natural number.  You input j = {j}.'.format(j=j))
        if not i < j:
            raise ValueError(
                'Index j must be greater than index i.  You input (i, j) = ({i}, {j}).'.format(i=i, j=j))
        seq = [0] * (max(i, j) + 1)
        seq[i] = 1
        seq[j] = -1
        seq = tuple(seq)
        return self._element_constructor_(seq)


class PieriOperatorAlgebra(ShiftingOperatorAlgebra):
    r""" The Pieri operator `u_i`.

    EXAMPLES::

        sage: u = PieriOperatorAlgebra()

    Act on catalan function::

        sage: cf = CatalanFunction([(0,2), (1,2)], [6, 6, 6])
        sage: u.i(2)(cf)
        CatalanFunction([(0,2), (1,2)], [6, 6, 5])

    Act on k-schur function::

        sage: base_ring = QQ['t']
        sage: Sym = SymmetricFunctions(base_ring)
        sage: t = base_ring.gen()
        sage: ks = Sym.kBoundedSubspace(4, t).kschur()
        sage: u.i(2)(ks[2, 2, 1])
        ks4[2, 2, 1] + t^2*ks4[3, 2] + t^3*ks4[4, 1]

    ..  SEEALSO::
        :class:`ShiftingOperatorAlgebra`
    """
    # TODO: verify by hand that above is really correct, or maybe a simpler example
    def __init__(self, base_ring=QQ['t'], prefix='u'):
        ShiftingOperatorAlgebra.__init__(self,
                                         base_ring=base_ring,
                                         prefix=prefix,
                                         basis_indices=ShiftingSequenceSpace())

    def i(self, i):
        r"""
        Return the Pieri operator `u_i`.

        Shorthand element constructor that allows you to create Pieri operators using the familiar `u_i` notation, with the exception that indices here are 0-based, not 1-based.

        EXAMPLES:

        Create the Pieri operator which lowers part 2 (indices are 0-based)::

            sage: u.i(2)
            u((0, 0, -1))

        ..  SEEALSO::
            :meth:`ShiftingOperatorAlgebra._element_constructor_`, :meth:`ShiftingOperatorAlgebra.__getitem__`
        """
        if not i in NonNegativeIntegerSemiring():
            raise ValueError(
                'i must be a natural number.  You input i = {i}.'.format(i=i))
        seq = [0] * (i + 1)
        seq[i] = -1
        seq = tuple(seq)
        return self._element_constructor_(seq)

    class Element(ShiftingOperatorAlgebra.Element):
        def __call__(self, operand):
            r"""
            Return the action of this raising operator element on the index ``operand``.

            EXAMPLES::

                sage: R = RaisingOperatorAlgebra()

                sage: R[1, -1]([1, 2, 3, 4])
                [([2, 1, 3, 4], 1)]

                sage: R[1, -1](([1, 2, 3, 4], [1, 1, 1]))
                ([([2, 1, 3, 4], 1)], [([2, 0, 1], 1)])

                sage: Sym = SymmetricFunctions(QQ['t'])
                sage: s = Sym.schur()
                sage: R[1, -1](s[4, 3, 2, 1])
                s[5, 2, 2, 1]
            """
            if _is_k_schur(operand):
                # convert to catalans
                kschur = operand.parent()
                base_ring = kschur.base_ring()
                cat_coeff_pairs = [
                    (CatalanFunctions().init_from_k_schur(
                        kschur(index), base_ring=base_ring), coeff)
                    for index, coeff in operand]
                # act
                new_cat_coeff_pairs = [(self.__call__(cat), coeff)
                                       for cat, coeff in cat_coeff_pairs]
                # convert back to kschur
                out_coeff_pairs = [(kschur(cat.eval()), coeff)
                                   for cat, coeff in cat_coeff_pairs]
                return sum(coeff * func for func, coeff in out_coeff_pairs)
            else:
                return ShiftingOperatorAlgebra.Element.__call__(self, operand)


class HallLittlewoodVertexOperator:
    r""" The Hall-Littlewood vertex operator.

    Garsia's version of Jing's Hall-Littlewood vertex operators.  These are defined in equations 4.2 and 4.3 of [cat]_ and appear visually as a bold capital H.

    INPUTS:

    - ``base_ring`` -- (default ``QQ['t']``) the base ring to build the :class:`SymmetricFunctions` upon.

    EXAMPLES::

        sage: H = HallLittlewoodVertexOperator
        sage: one = SymmetricFunctions(QQ['t']).hall_littlewood().Qp().one()
        sage: H([4, 1, 3])(one) == H(4)(H(1)(H(3)(one)))
        True

    ..  SEEALSO::
        :meth:`sage.combinat.sf.new_kschur.KBoundedSubspaceBases.ElementMethods.hl_creation_operator`
    """

    def __init__(self, composition, base_ring=QQ['t'], t=None):
        if composition in NonNegativeIntegerSemiring():
            self.composition = Composition([composition])
        elif _is_sequence(composition):
            self.composition = Composition(composition)
        else:
            raise ValueError('Bad composition.')
        self.base_ring = base_ring
        self.sym = SymmetricFunctions(self.base_ring)
        self._t = t

    def __repr__(self):
        r""" Return a human-friendly string representation of this Hall-Littlewood vertex operator.

        This string also serves as an example of what somebody might type to create this Hall-Littlewood vertex operator in sage.

        EXAMPLES::

            sage: H = HallLittlewoodVertexOperator
            sage: H([4, 1, 3])
            HallLittlewoodVertexOperator([4, 1, 3])
        """
        return '{}({})'.format(self.__class__.__name__, self.composition)

    def _hh(self, k):
        r""" Internal helper function.

        Return the homogeneous function indexed by the integer ``k``.

        EXAMPLES::

            sage: op = HallLittlewoodVertexOperator([4, 1, 3])
            sage: op._hh(-2)
            0
            sage: op._hh(-1)
            0
            sage: op._hh(0)
            s[]
            sage: op._hh(1)
            s[1]
            sage: op._hh(2)
            s[2]
        """
        # homogeneous indexed by an integer k (positive or negative)
        # if k is less than 0, result is 0
        # if k == 0, result is s([])
        # if k > 0, then the result is s([k])
        sym = self.sym
        s = sym.s()
        if k == 0:
            return s.one()
        elif k < 0:
            # 0, but as a sym func
            return 0 * s.one()
        elif k > 0:
            return s([k])
        else:
            raise ValueError

    def _skewbyeeq(self, k, f):
        r""" Internal helper function.

        Given integer ``k`` and function ``f``, return the function `f` skewed by e[k](1-t).
        """
        # skew by e[k](1-t)
        sym = self.sym
        e = sym.e()
        if self._t is None:
            t = self.base_ring.gen()
        else:
            t = self._t
        if k == 0:
            return f
        elif k < 0:
            # 0, but as a sym func
            return 0 * f
        elif k > 0:
            return f.skew_by(e([k]).theta_qt(t, 0))
        else:
            raise ValueError

    def _op(self, m, f):
        r""" Internal helper function.

        This method is Jing's Hall-Littlewood creation operator.

        ..  SEEALSO::
            :meth:`sage.combinat.sf.new_kschur.KBoundedSubspaceBases.ElementMethods.hl_creation_operator`
        """
        return sum((-1)**k * self._hh(m+k) * self._skewbyeeq(k, f) for k in range(f.degree() + 1))

    def __call__(self, input_):
        r""" Return the action of this Hall-Littlewood operator on ``input_``.

        INPUTS:

        - ``input_`` -- A symmetric function.

        EXAMPLES:

        We typically act on the identity, so let's retrieve the identity::

            sage: Sym = SymmetricFunctions(QQ['t'])
            sage: hl = Sym.hall_littlewood().Qp()
            sage: one = hl.one()

        Let's look at the action of `H_{4, 1, 3}` on the identity::

            sage: H = HallLittlewoodVertexOperator
            sage: H([4, 1, 3])(one)
            (t-1)*HLQp[4, 2, 2] + t*HLQp[4, 3, 1]

        ..  SEEALSO::
            :meth:`sage.combinat.sf.new_kschur.KBoundedSubspaceBases.ElementMethods.hl_creation_operator`
        """
        gamma = self.composition
        sym = self.sym
        if self._t is None:
            HLQp = sym.hall_littlewood().Qp()
        else:
            HLQp = sym.hall_littlewood(t=self._t).Qp()
        # iterate
        for part in reversed(gamma):
            input_ = self._op(part, input_)
        return HLQp(input_)


def compositional_hall_littlewood_Qp(gamma, base_ring=QQ['t'], t=None):
    r""" Given gamma, returns the compositional Hall-Littlewood polynomial `H_{\gamma}(\mathbf{x}; t)` in the Q' basis, as defined in [cat]_ section 4.4.

    If the composition gamma is a partition, this is just the Hall-Littlewood Q' polynomial.

    EXAMPLES::

        sage: hl = SymmetricFunctions(QQ['t']).hall_littlewood().Qp()
        sage: compositional_hall_littlewood_Qp([3, 3, 2]) == hl[3, 3, 2]
        True

    ..  SEEALSO::
        :meth:`sage.combinat.sf.hall_littlewood.HallLittlewood.Qp`
    """
    sym = SymmetricFunctions(base_ring)
    if t is None:
        HLQp = sym.hall_littlewood().Qp()
    else:
        HLQp = sym.hall_littlewood(t=t).Qp()
    if is_weakly_decreasing(gamma) and all(term > 0 for term in gamma[:-1]) and ((not gamma) or gamma[-1] >=0):
        # this is MUCH faster than the HallLittlewoodVertexOperator for partitions of length 5ish
        gamma = Partition(gamma)
        return HLQp(gamma)
    else:
        H = HallLittlewoodVertexOperator
        return H(gamma, base_ring=base_ring, t=t)(HLQp.one())

def raising_roots_operator(roots, base_ring=QQ['t'], t=1):
    r""" Return the operator `\prod_{(i,j) \in roots} (1 - tR_{ij})`.

    Given a list of roots `roots = \Phi` (often a root ideal), and optionally a variable `t`, return the operator

    ..  MATH::

        \prod_{(i,j) \in \Phi} (1 - tR_{ij}).

    If you input an integer for roots (e.g. ``roots = 3``), it will use the biggest possible root ideal in the `n` x `n` grid (the '`n`-th staircase root ideal').

    EXAMPLES:

    Consider the root ideal [(0, 1)] in the 2 x 2 grid.  This root ideal gives the raising roots operator `1 - R_{01}`::

        sage: R = RaisingOperatorAlgebra()
        sage: 1 - R.ij(0, 1)
        R() - R(1, -1)

        sage: raising_roots_operator([(0, 1)])
        R() - R(1, -1)

        sage: raising_roots_operator(2)
        R() - R(1, -1)

    We can pass in `t` in addition to get the raising roots operator `1 - t R_{01}`.  We can also choose a base ring to work over::

        sage: K = FractionField(ZZ['t'])
        sage: t = K.gen()
        sage: R = RaisingOperatorAlgebra(base_ring=K)

        sage: 1 - t * R.ij(0, 1)
        R() - t*R(1, -1)

        sage: raising_roots_operator([(0, 1)], base_ring=K, t=t)
        R() - t*R(1, -1)

        sage: raising_roots_operator(2, base_ring=K, t=t)
        R() - t*R(1, -1)

    ..  SEEALSO::
        :meth:`qt_raising_roots_operator`, :class:`CatalanFunction`
    """
    R = RaisingOperatorAlgebra(base_ring=base_ring)
    if roots in NonNegativeIntegerSemiring():
        roots = RootIdeals().init_staircase(roots)
    op = prod([1 - t*R.ij(i, j) for (i, j) in roots], R.one())
    return op


def qt_raising_roots_operator(roots, base_ring=QQ['t', 'q'], t=None, q=None):
    r""" Return the operator `\prod_{ij \in \Phi} (1 - tR_{ij}) \prod_{ij \in roots} (1 - qR_{ij})`.

    The q-t analogue of :meth:`raising_roots_operator`, defined by

    ..  MATH::
        \prod_{ij \in \Phi} (1 - tR_{ij}) \prod_{ij \in \Phi} (1 - qR_{ij}).

    EXAMPLES:

    Consider the root ideal [(0, 1)] in the 2 x 2 grid.  This root ideal gives the "qt" raising roots operator `(1 - t R_{01})(1 - q R_{01})` which equals `1 + (-t-q) R_{01} + t q R_{01}^2`::

        sage: K = QQ['t', 'q']
        sage: (t, q) = K.gens()
        sage: R = RaisingOperatorAlgebra(base_ring=K)

        sage: (1 - t*R.ij(0, 1)) * (1 - q*R.ij(0, 1))
        R() + (-t-q)*R(1, -1) + t*q*R(2, -2)

        sage: qt_raising_roots_operator([(0, 1)], base_ring=K, t=t, q=q)
        R() + (-t-q)*R(1, -1) + t*q*R(2, -2)

    Since ``QQ['t', 'q']`` happens to be the default base ring, and its generators are ``t`` and ``q``, we did not actually need to pass in the base ring, `t`, nor `q` in the above example::

        sage: qt_raising_roots_operator([(0, 1)])
        R() + (-t-q)*R(1, -1) + t*q*R(2, -2)

    ..  SEEALSO::
        :meth:`raising_roots_operator`
    """
    if roots in NonNegativeIntegerSemiring():
        roots = RootIdeals().init_staircase(roots)
    if t is None:
        t = base_ring.gens()[0]
    if q is None:
        q = base_ring.gens()[1]
    op1 = raising_roots_operator(roots, base_ring=base_ring, t=q)
    op2 = raising_roots_operator(roots, base_ring=base_ring, t=t)
    return op2 * op1


class CatalanFunction:
    r""" A catalan function `H(\Psi; \gamma)` as discussed in [cat]_ section 4.

    By definition,

    ..  MATH::

        H(\Psi; \gamma) := \prod_{ij \in \Psi} (1 - R_{ij})^{-1} s_\gamma = \prod_{ij \in \Delta^+ \smallsetminus \Psi} (1 - R_{ij}) H_\gamma

    where `s_\gamma` is a schur function and `H_\gamma` is a Hall-Littlewood Q' function.

    The parameters for initializing a catalan function are exactly the same as in :meth:`CatalanFunctions.init_from_indexed_root_ideal`.

    EXAMPLES::

        sage: CatalanFunction([(0,2), (1,2)], [6, 6, 5])
        H([(0, 2), (1, 2)]; [6, 6, 5])

    ..  SEEALSO::
        There are more ways to initialize a catalan function, and the methods for doing so are found in :class:`CatalanFunctions`.
    """
    BASE_RING_DEFAULT = QQ['t']
    PREFIX_DEFAULT = 'H'

    def __init__(self, roots, index, base_ring=None, prefix=None):
        assert _is_sequence(index)
        self.index = index
        assert root_ideal._is_roots(roots)
        self.roots = RootIdeal(roots, len(self.index))
        self.base_ring = base_ring if base_ring is not None else self.BASE_RING_DEFAULT
        self.prefix = prefix if prefix is not None else self.PREFIX_DEFAULT

    def __eq__(self, other):
        r""" Return whether this catalan function is equal to the catalan function ``other``.

        EXAMPLES::

            sage: cf1 = CatalanFunction([(0,2), (1,2)], [6, 6, 5])
            sage: cf2 = CatalanFunction({(1,2), (0,2)}, Partition([6, 6, 5]))
            sage: cf1 == cf2
            True
        """
        # TODO: account for the fact that DIFFERENT root/index pairs could actually give the SAME catalan function!!
        return self.roots == other.roots and self.index == other.index and self.base_ring == other.base_ring

    def __repr__(self):
        r""" Return a human-friendly string representation of this catalan function.

        EXAMPLES::

            sage: cf = CatalanFunction([(0,2), (1,2)], [6, 6, 5])
            sage: cf
            H([(0, 2), (1, 2)]; [6, 6, 5])
        """
        return '{}({}; {})'.format(self.prefix, self.roots, self.index)

    def _get_indices_for_index_operator(self):
        r""" Internal helper function.

        This method is defined so that an element of any :class:`ShiftingOperatorAlgebra` (such as a raising operator or a pieri operator) can act on this catalan function.

        EXAMPLES::

            sage: cf = CatalanFunction([(0,2), (1,2)], [6, 6, 5])
            sage: cf._get_indices_for_index_operator()
            [6, 6, 5]

        This is what allows you to do something like::

            sage: cf = CatalanFunction([(0,2), (1,2)], [6, 6, 5])
            sage: R = RaisingOperatorAlgebra()
            sage: R.ij(0, 1)(cf)
            H([(0, 2), (1, 2)]; [7, 5, 5])

        ..  SEEALSO::
            :meth:`_new_object_for_index_operator`
        """
        return self.index

    def _new_object_for_index_operator(self, new_index):
        r""" Internal helper function.

        This method is defined so that an element of any :class:`ShiftingOperatorAlgebra` (such as a raising operator or a pieri operator) can act on this catalan function.

        Returns a brand new :class:`CatalanFunction`.

        EXAMPLES::

            sage: cf = CatalanFunction([(0,2), (1,2)], [6, 6, 5])
            sage: cf._new_object_for_index_operator([7, 5, 5])
            H([(0, 2), (1, 2)]; [7, 5, 5])

        This is what allows you to do something like::

            sage: cf = CatalanFunction([(0,2), (1,2)], [6, 6, 5])
            sage: R = RaisingOperatorAlgebra()
            sage: R.ij(0, 1)(cf)
            H([(0, 2), (1, 2)]; [7, 5, 5])

        ..  SEEALSO::
            :meth:`_get_indices_for_index_operator`
        """
        new_obj = self.__class__(self.roots, new_index,
                                 base_ring=self.base_ring)
        return new_obj

    def eval(self, t=None):
        r""" Return the catalan function in terms of the Hall-Littlewood Q' basis.

        EXAMPLES::

            sage: delta_plus = partition_to_root_ideal([2, 1], n=3)
            sage: cf = CatalanFunction(delta_plus, [3, 1, 1])
            sage: cf.eval()
            HLQp[3, 1, 1]

            sage: s = SymmetricFunctions(QQ['t']).schur()
            sage: cf = CatalanFunction([], [4, 1])
            sage: cf.eval()
            HLQp[4, 1] - t*HLQp[5]
            sage: s(cf.eval())
            s[4, 1]

        This function can also specialize `t` to an arbitrary integer::

            sage: cf.eval(t=1)
            HLQp[4, 1] - HLQp[5]
            sage: cf.eval(t=-1)
            HLQp[4, 1] + HLQp[5]
            sage: s(cf.eval(t=-1))
            s[4, 1]

        ..  SEEALSO::
            :meth:`expand`
        """
        # setup
        if t is None:
            hl = SymmetricFunctions(self.base_ring).hall_littlewood().Qp()
            t = self.base_ring.gen()
        else:
            hl = SymmetricFunctions(self.base_ring).hall_littlewood(t=t).Qp()
        # formula
        roots_complement = self.roots.complement()
        op = raising_roots_operator(
            roots_complement, base_ring=self.base_ring, t=t)
        hl_poly = hl(self.index)
        cat_func = op(hl_poly)
        return cat_func

    def expand(self, *args, **kwargs):
        r"""
        Expand the catalan function as a symmetric polynomial in ``n`` variables.

        INPUT:

        - ``n`` -- a nonnegative integer

        - ``alphabet`` -- (default: ``'x'``) a variable for the expansion

        OUTPUT:

        A monomial expansion of ``self`` in the `n` variables
        labelled ``x0``, ``x1``, ..., ``x{n-1}`` (or just ``x``
        if `n = 1`), where ``x`` is ``alphabet``.

        EXAMPLES::

            sage: cf = CatalanFunction([], [4, 1])
            sage: cf.expand(1)
            0
            sage: cf.expand(2)
            x0^4*x1 + x0^3*x1^2 + x0^2*x1^3 + x0*x1^4
            sage: cf.expand(2, alphabet='y')
            y0^4*y1 + y0^3*y1^2 + y0^2*y1^3 + y0*y1^4

        ..  SEEALSO::
            :meth:`eval`
        """
        return self.eval().expand(*args, **kwargs)

    def _latex_(self, color='red'):
        r"""
        Return LaTeX code to draw a LaTeX representation of a root ideal 
        encoding  ``self``.

        EXAMPLES::

            sage: cf = CatalanFunction([(0,2)],[3,2,1])
            sage: latex(cf) # indirect doctest
            \begin{ytableau}
              3 & {} & *(red) \\ 
              {} & 2 & {} \\ 
              {} & {} & 1 
            \end{ytableau}
        """
        return self.roots._latex_(color=color,index=self.index)
    
    # def _latex_(self):
    #     r"""
    #     Return LaTeX code to draw a LaTeX representation of a root ideal 
    #     encoding  ``self``.

    #     EXAMPLES::

    #         sage: cf = CatalanFunction([(0,1)], [1,2])
    #         sage: latex(cf) # indirect doctest
    #         \begin{tikzpicture}[every node/.style={minimum size=.5cm-\pgflinewidth, outer sep=0pt}] 
    #         \draw[step=0.5cm,color=black] (0,0) grid (1.0,-1.0); 
    #           \node[fill=red] at (0.75,-0.25) {}; 
    #           \node at (0.25,-0.25) {1}; 
    #           \node at (0.75,-0.75) {2}; 
    #         \end{tikzpicture}
    #     """
    #     ideal_tikz = latex(self.roots).split("\n")
    #     index_nodes = ["  \\node at (" + str(i/2.0+0.25) + "," + str(-1*(i/2.0+0.25)) + ") {" + str(self.index[i]) + "};" for i in range(len(self.index))]
    #     all_lines = ideal_tikz[:-1] + index_nodes + [ideal_tikz[-1]]
    #     result = reduce(lambda a,b: a+"\n"+b, all_lines)
    #     return result
        
class CatalanFunctions:
    r""" The family of catalan functions, as discussed in [cat]_ section 4.

    Use this class as a factory to initialize a :class:`CatalanFunction` object with any valid identifying data.  See the ``init_from...`` methods below for all possible ways to create a catalan function.
    """

    def init_from_indexed_root_ideal(self, roots, index, base_ring=None, prefix=None):
        r"""
        INPUTS:

        - ``roots`` -- iterable of roots `\Phi` (typically a root ideal)

        - ``index`` -- composition `\gamma` that indexes the root ideal and appears in `s_\gamma` and `H_\gamma` below

        OPTIONAL INPUTS:

        - ``base_ring`` -- (default ``QQ['t']``) the ring over which to build the `h_\gamma(x; \alpha)`'s

        OUTPUTS:

        The catalan function

        ..  MATH::

            H(\Phi; \gamma) := \prod_{ij \in \Phi} (1 - R_{ij})^{-1} s_\gamma = \prod_{ij \in \Delta^+ \smallsetminus \Phi} (1 - R_{ij}) H_\gamma

        where `s_\gamma` is a schur function and `H_\gamma` is a Hall-Littlewood Q' function.

        EXAMPLES::

            sage: CatalanFunctions().init_from_indexed_root_ideal([(0,2), (1,2)], [6, 6, 5])
            H([(0, 2), (1, 2)]; [6, 6, 5])
        """
        return CatalanFunction(roots, index, base_ring, prefix)

    def init_from_skew_partition(self, sp, base_ring=None, prefix=None, algorithm='removable roots'):
        r""" Given a skew partition ``sp``, return the catalan function `H(\Phi^+(sp); rs(sp))`.

        Given a linked :class:`SkewPartition` ``sp``, return the :class:`CatalanFunction` `H(\Phi^+(sp); rs(sp))`.

        EXAMPLES::

            sage: CFS = CatalanFunctions()
            sage: sp = SkewPartition([[2, 1, 1], [1]])
            sage: CFS.init_from_skew_partition(sp)
            H([(0, 1), (0, 2)]; [1, 1, 1])

        ..  SEEALSO::
            :meth:`init_from_row_and_column_lengths`
        """
        ri = RootIdeals().init_from_skew_partition(sp, type='max', algorithm=algorithm)
        rs = sp.row_lengths()
        return self.init_from_indexed_root_ideal(ri, rs, base_ring, prefix)

    def init_from_row_and_column_lengths(self, row_lengths, column_lengths, base_ring=None, prefix=None, algorithm='removable roots'):
        r""" Determine the skew partition `D` with row-shape ``row_lengths`` and column-shape ``column_lengths``, and return the catalan function `H(\Phi^+(D); \text{row_lengths})`.

        EXAMPLES::

            sage: CFS = CatalanFunctions()
            sage: CFS.init_from_row_and_column_lengths([1, 1, 1], [2, 1])
            H([(0, 1), (0, 2)]; [1, 1, 1])

        ..  SEEALSO::
            :meth:`init_from_skew_partition`
        """
        sp = SkewPartitions().from_row_and_column_length(row_lengths, column_lengths)
        return self.init_from_skew_partition(sp, base_ring=base_ring, prefix=prefix, algorithm=algorithm)

    def init_from_k_shape(self, p, k, base_ring=None, prefix=None, algorithm='removable roots'):
        r""" Given ``k`` and a `k`-shape ``p``, return the catalan function `H(\Phi^+(p); rs(p))`.

        EXAMPLES::

            sage: CFS = CatalanFunctions()
            sage: CFS.init_from_k_shape([4, 3, 2, 1], 1)
            H([]; [4, 3, 2, 1])
        """
        assert is_k_shape(p, k)
        sp = SkewPartition([p, []])
        return self.init_from_skew_partition(sp, base_ring=base_ring, prefix=prefix, algorithm=algorithm)

    def init_from_k_schur(self, func, base_ring=None, prefix=None):
        r""" Given a `k`-schur function ``func`` `= s^k_\lambda(x;t)`, initialize the catalan function `H(\Delta^k(\lambda); \lambda)`.

        Mathematically, these two functions are equal.  The usefulness of this method is that you input a ``sage.combinat.sf.new_kschur.kSchur_with_category`` object and you obtain a :class:`CatalanFunction` object.

        EXAMPLES::

            sage: base_ring = QQ['t']
            sage: Sym = SymmetricFunctions(base_ring)
            sage: t = base_ring.gen()
            sage: ks = Sym.kBoundedSubspace(4, t).kschur()
            sage: func = ks[2, 1, 1]
            sage: CatalanFunction(func, base_ring=base_ring)
            H([]; [2, 1, 1])
        """
        # TODO: make sure [] above is correct.  go by hand.
        # check inputs
        assert _is_k_schur(func)
        assert len(func.support()) == 1
        # gather roots
        index = func.support()[0]
        k = func.parent().k
        roots = partition_to_k_schur_root_ideal(index, k)
        # initialize from roots and index
        return self.init_from_indexed_root_ideal(roots, index, base_ring, prefix)

    def init_parabolic_from_composition_and_index(self, composition, index, base_ring=None, prefix=None):
        r""" Given a composition `\eta` of positive integers and an index `\gamma`, return the parabolic catalan function `H(\Delta(\eta), \gamma)`, where

        ..  math::

            \Delta(\eta) := \{ \alpha \in \Delta^+_{|\eta|} \:\text{above the block diagonal with block sizes}\: \eta_1, \ldots, \eta_r\}

        as in [cat]_ just below Conjecture 3.3.

        EXAMPLES::

            sage: CatalanFunctions().init_parabolic_from_composition_and_index([1, 3, 2], [1, 2, 3, 4, 5, 6])
            H([(0, 1), (0, 2), (0, 3), (0, 4), (0, 5), (1, 4), (1, 5), (2, 4), (2, 5), (3, 4), (3, 5)]; [1, 2, 3, 4, 5, 6])
        """
        ri = RootIdeals().init_parabolic_from_composition(composition)
        return self.init_from_indexed_root_ideal(ri, index, base_ring, prefix)


##############
def k_plus_one_core_to_k_schur_function(p, k, base_ring=QQ['t']):
    r""" Given a `k+1`-core ``p``, the natural number ``k`` itself, and optionally a ``base_ring``, return the corresponding `k`-schur function. """
    # TODO: compare the performance of this function to existing k-schur function.
    p = Core(p, k+1)
    return k_shape_to_catalan_function(p, k, base_ring)


class InfiniteDimensionalFreeAlgebra(CombinatorialFreeModule):
    r"""
    By default, the algebra generated by ``x[0], x[1], x[2], ...`` over the integers.

    To change the index set of the generators, use ``index_set=`` (default ``NN``).  To overhaul the set of generators entirely (not recommended), use ``basis_indices=``.

    To change the ring that the algebra works over, use ``base_ring=`` (default ``ZZ``).

    To change the prefix of the generators, use ``prefix=`` (default ``'x'``).

    ..  SEEALSO::
        sage.algebras.free_algebra, https://doc.sagemath.org/html/en/reference/monoids/sage/monoids/free_abelian_monoid.html, sage.monoids.free_monoid, sage.monoids.free_abelian_monoid
    """

    def __init__(self,
                 base_ring=IntegerRing(),
                 prefix='x',
                 basis_indices=None,
                 index_set=NonNegativeIntegerSemiring()):
        self._base_ring = base_ring
        self._index_set = index_set
        self._basis_monoid = FreeMonoid(
            index_set=self._index_set, commutative=True, prefix=prefix) if basis_indices is None else basis_indices
        # category
        category = Algebras(self._base_ring.category()
                            ).WithBasis().Commutative()
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
        r""" Return whether the infinite dimensional free algebra is a prime field (it never is!).

        Note that this method exists for internal compatability.

        EXAMPLES::

            sage: A = InfiniteDimensionalFreeAlgebra()
            sage: A.is_prime_field()
            False
        """
        return False

    def _element_constructor_(self, monoid_el):
        r""" Return the element whose index is ``monoid_el``.

        Note: This method is only for basis elements.  This method for internal use.

        EXAMPLES::

        ..  SEEALSO::
            :meth:`__getitem__`
        """
        assert monoid_el in self._basis_monoid
        return self.basis()[monoid_el]

    def __getitem__(self, user_input):
        r""" Given an integer ``user_input``, return the corresponding basis element.

        EXAMPLES::

            sage: A = InfiniteDimensionalFreeAlgebra()
            sage: A[4]
            x[4]

        ..  SEEALSO::
            :meth:`_element_constructor_`
        """
        # USER front entrance to creating elements "x[4]"
        assert user_input in IntegerRing()
        monoid_el = self._basis_monoid.gen(user_input)
        return self.basis()[monoid_el]

    @cached_method
    def one_basis(self):
        r""" Return the index of the identity element.

        EXAMPLES::

            sage: A.one_basis()
            1
        """
        return self._basis_monoid.one()

    def product_on_basis(self, monoid_el1, monoid_el2):
        r""" Given indices ``monoid_el1`` and ``monoid_el2``, return the product of the basis elements indexed by those indices.



        """
        monoid_el_product = monoid_el1 * monoid_el2
        return self._element_constructor_(monoid_el_product)

    def _repr_(self):
        r""" Return a human-friendly string representation of this infinite-dimensional free algebra.

        EXAMPLES::

            sage: F = InfiniteDimensionalFreeAlgebra()
            sage: F
            InfiniteDimensionalFreeAlgebra_with_category with generators indexed by Non negative integer semiring, over Integer Ring

            sage: F = InfiniteDimensionalFreeAlgebra(base_ring=QQ, index_set=ZZ)
            sage: F
            InfiniteDimensionalFreeAlgebra_with_category with generators indexed by Integer Ring, over Rational Field
        """
        return "{class_name} with generators indexed by {index_set}, over {base_ring}".format(class_name=self.__class__.__name__, base_ring=self._base_ring, index_set=self._index_set)


DoubleRing = InfiniteDimensionalFreeAlgebra(prefix='a', index_set=IntegerRing())
r""" ``DoubleRing`` is the ring `\Lambda(a)` found in [Fun]_ section 3.

EXAMPLES::

    sage: DoubleRing
    InfiniteDimensionalFreeAlgebra_with_category with generators indexed by Integer Ring, over Integer Ring

..  SEEALSO::
    :class:`InfiniteDimensionalFreeAlgebra`
"""


def dual_k_theoretic_homogeneous(k, r, base_ring=QQ):
    r""" The dual K-theoretic h, often denoted Kh, is defined for any integer `k` by the formula `h_k(x, r) = \sum_{i=0}^{k} \binom{r + i - 1}{i} h_{k - i}(x)` in [LN]_ p.88 top-right.

    If ``k`` and ``r`` are compositions, then it is recursively defined as `h_k(x, r) := \prod_j h_{k_j}(x, r_j)`.

    EXAMPLES::

        sage: dual_k_theoretic_homogeneous(0, 0)
        1

        sage: dual_k_theoretic_homogeneous(1, 2, base_ring=QQ['t'])
        h[1] + 2

        sage: dual_k_theoretic_homogeneous([2, 1], [1, 1])
        h[1]**2 + h[1]*h[2] + 2*h[1] + h[2] + 1

    ..  SEEALSO::
        :meth:`dual_k_catalan_function`
    """
    if _is_sequence(k):
        # pad with 0's
        max_len = max(len(k), len(r))
        k = list(k) + [0] * (max_len - len(k))
        r = list(r) + [0] * (max_len - len(r))
        # multiply
        h_list = [dual_k_theoretic_homogeneous(
            k_el, r_el, base_ring) for k_el, r_el in zip(k, r)]
        return reduce(operator.mul, h_list)
    else:
        assert k >= 0
        h = SymmetricFunctions(base_ring).h()
        return sum(binomial(r + i - 1, i) * h[k - i] for i in range(k + 1))


def dual_k_catalan_function(roots, index, index2, base_ring=QQ):
    r"""
    INPUTS:

    - ``roots`` -- iterable of roots `\Phi` (typically a :class:`RootIdeal`)

    - ``index`` -- composition `\gamma` that indexes `h_\gamma(x; \alpha)`

    - ``index2`` -- composition `\alpha` used in `h_\gamma(x; \alpha)`

    OPTIONAL INPUTS:

    - ``base_ring`` -- (default ``QQ``) the ring over which to build the `h_\gamma(x; \alpha)`'s

    OUTPUTS:

    The 'dual k' catalan function

    ..  MATH::

        \prod_{ij \in \Delta^+ \smallsetminus \Phi} (1 - R_{ij}) h_\gamma(x; \alpha).

    ..  SEEALSO::
        :meth:`dual_k_theoretic_homogeneous`, :meth:`dual_grothendieck_function`
    """
    # setup
    Kh = dual_k_theoretic_homogeneous(index, index2, base_ring=base_ring)
    # formula
    roots = RootIdeal(roots, n=len(index))
    roots_complement = roots.complement()
    op = raising_roots_operator(roots_complement, base_ring=base_ring, t=1)
    cat_func = op(Kh)
    return cat_func


def dual_grothendieck_function(composition, base_ring=QQ):
    r""" Given a composition `composition = \lambda`, return the dual Grothendieck function defined by `g_\lambda(x) = \text{det}(h_{\lambda_i + j - i}(x, i - 1))` in [LN]_ p.88 equation (4).

    Equivalently, the dual Grothendieck function is `g_\lambda(x) = \prod_{ij \in \Delta^+} (1 - R_{ij}) h_\lambda(x; (0, 1, \ldots, n-1))`.

    EXAMPLES::

        sage: h = SymmetricFunctions(QQ).h()
        sage: dual_grothendieck_function([2, 1])
        h[1]*h[2] + h[2] - h[3]

    ..  SEEALSO::
        :meth:`dual_k_catalan_function`
    """
    roots = []  # because dual_k_catalan_function will take the complement
    n = len(composition)
    reversed_staircase_ptn = list(reversed(staircase_shape(n)))
    return dual_k_catalan_function(roots, composition, reversed_staircase_ptn, base_ring=base_ring)


def double_homogeneous_building_block(p, n):
    r""" The double complete homogeneous symmetric polynomial "building block" `h_p(x || a)`.

    Defined as

    ..  MATH::
        h_p(x_1, \ldots, x_n \,||\, a) \,= \sum_{n \geq i_1 \geq \ldots \geq i_p \geq 1} (x_{i_1} - a_{i_1})(x_{i_2} - a_{i_2 - 1}) \cdots (x_{i_p} - a_{i_p - p + 1})

    in [Fun]_ section 3 between equation (6) and (7).  Note that our indices are 0-based.

    ..  SEEALSO::
        :class:`DoubleHomogeneous`, :meth:`double_homogeneous_building_block_shifted`, :meth:`shift`
    """
    a = DoubleRing
    sym = SymmetricFunctions(DoubleRing)
    s = sym.s()
    one = s.one()
    one_poly = one.expand(n)
    x = one_poly.parent().gens()
    ptns = Partitions(max_part=n, length=p)
    total_sum = 0
    for ptn in ptns:
        summand = prod([(x[ptn[b]] - a[ptn[b] - b]) for b in range(p)], a.one())
        total_sum += summand
    return total_sum


def shift(element):
    r""" The function `\tau` which acts on any element of `\Lambda(a)` (``DoubleRing``) by sending each element `a_i` to `a_{i+1}` for all `i`.  It can be found in [Fun]_ p.8 between equations (6) and (7).

    ..  SEEALSO::
        :meth:`double_homogeneous_building_block_shifted`, :meth:`double_homogeneous_building_block`
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

    ..  SEEALSO::
        :class:`DoubleHomogeneous`, :meth:`double_homogeneous_building_block`, :meth:`shift`
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
    INPUTS:

    ``mu1`` -- composition

    ``mu2`` -- composition

    ``n`` -- number of `x` variables

    EXAMPLES:

    Create the double homogeneous `h^{(4)}_{\mu, \beta}`::

        sage: DoubleHomogeneous(mu, beta, 4)

    Create the double homogeneous shifted building block `h_{r, s}` in 4 variables::

        sage: r = 5
        sage: s = 2
        sage: DoubleHomogeneous([r], [s], 4)

    Create the double homogeneous symmetric building block `h_p(x \,||\, a)` in 4 variables::

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
        r""" Internal helper function.

        This method is defined so that an element of any :class:`ShiftingOperatorAlgebra` (such as a raising operator) can act on this double homogeneous function.
        """
        pair = (self.index1, self.index2)
        return pair

    def _new_object_for_index_operator(self, indices):
        r""" Internal helper function.

        This method is defined so that an element of any :class:`ShiftingOperatorAlgebra` (such as a raising operator) can act on this double homogeneous function.
        """
        (index1, index2) = indices
        new_obj = self.__class__(index1, index2, self.n)
        return new_obj

    def eval(self, n):
        r""" Given the number of variables ``n``, return ``self`` expanded in terms of the shifted double homogeneous building blocks `h_{r, s}`. """
        (mu1, mu2) = (self.index1, self.index2)
        # pad with 0's
        max_len = max(len(mu1), len(mu2))
        mu1 = list(mu1) + [0] * (max_len - len(mu1))
        mu2 = list(mu2) + [0] * (max_len - len(mu2))
        # compute
        h_list = [double_homogeneous_shifted(
            mu1[i], mu2[i], n) for i in range(max_len)]
        hp = prod(h_list)
        return hp

    def expand(self):
        raise NotImplemented


def double_schur(index, n):
    r""" Given a composition ``index`` `= \lambda` and the number of variables `n`, return the double Schur function defined by

    ..  MATH::
        s_\lambda(x_1, \ldots, x_n \,||\, a) = \text{det}\left(h^{(n)}_{\lambda_i + i - j, j - 1}\right)

    or equivalently, defined by

    ..  MATH::
        s_\lambda(x_1, \ldots, x_n \,||\, a) = \prod_{ij \in \Delta^+(l)} (1 - R_{ij}) h^{(n)}_{\lambda, (0, \ldots, l-1)}

    where `l` is the length of `\lambda`, in [Fun]_ p.9 equation (9).
    """
    l = len(index)
    op = raising_roots_operator(l, base_ring=ZZ, t=1)
    rho = list(reversed(staircase_shape(l)))
    h_index = DoubleHomogeneous(index, rho, n)
    return op(h_index)


def double_catalan_function(roots, index, n):
    r""" Given some ``roots`` (typically a :class:`RootIdeal`), an ``index``, and the something something of the double homemgeneous function ``n``, return the corresponding double catalan function.
    """
    # TODO: test this
    # setup
    l = len(index)
    rho = list(reversed(staircase_shape(l)))
    h_index = DoubleHomogeneous(index, rho, n)
    # formula
    roots_complement = root_ideal.complement(roots, l)
    # TODO: see what to pass in for base_ring in below line.
    op = raising_roots_operator(roots_complement, t=1)
    cat_func = op(h_index)
    return cat_func


def substitute(f, t=None, q=None):
    r""" Given a symmetric function ``f``, plug the inputted ``t`` and ``q`` values into it and return the resulting function.

    ..  SEEALSO::
        :meth:`SymmetricFunctionAlgebra_schur_with_category.element_class.plethysm`, :meth:`ungraded`
    """
    # a t value of 'None' will leave t as-is.
    basis = f.parent()
    base_ring = f.base_ring()
    s = SymmetricFunctions(base_ring)
    f_s = s(f)
    coeffs = f_s.coefficients()
    # Necessary because otherwise coeffs and monomials don't line up
    monomials = sorted(f_s.monomials())
    specialized_coeffs = [coeff.substitute(t=t, q=q) for coeff in coeffs]
    combine = zip(specialized_coeffs, monomials)
    ungraded_f_s = sum(coeff * monom for (coeff, monom) in combine)
    ungraded_f = basis(ungraded_f_s)
    return ungraded_f


def ungraded(f):
    r""" Given a symmetric function ``f``, return the result of plugging in `t = 1`.

    ..  SEEALSO::
        :meth:`substitute`
    """
    return substitute(f, t=1)
