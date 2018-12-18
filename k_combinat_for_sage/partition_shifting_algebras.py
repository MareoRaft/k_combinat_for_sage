# -*- coding: utf-8 -*-
r"""
This module contains families of operators that act on partitions or. more generally, integer sequences. In particular, this includes Young's raising operators, which act on partitions by moving a cell from a smaller row to a larger row, thus resulting in a partition that is higher in dominance order. 

AUTHORS:

- Matthew Lancellotti, George H. Seelinger (2018): Initial version

REFERENCES:

.. [Fun] `Raising operators and the Littlewood-Richardson polynomials <https://arxiv.org/pdf/1203.4729.pdf>`_.  Fun, Alex.
.. [LN] `Finite sum Cauchy identity for dual Grothendieck polynomials <https://projecteuclid.org/download/pdf_1/euclid.pja/1407415930>`_.
"""

#*****************************************************************************
#  Copyright (C) 2018 Matthew Lancellotti <mvlancellotti@gmail.com>
#                     George H. Seelinger <ghseeli@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.algebras.group_algebra import GroupAlgebra_class
from sage.categories.algebras import Algebras
from sage.categories.groups import Groups
#from sage.combinat.sf.sf import SymmetricFunctions
from sage.combinat.composition import Composition
from sage.combinat.free_module import CombinatorialFreeModule
from sage.combinat.partition import Partitions, Partition
from sage.combinat.sf.sf import SymmetricFunctions
from sage.combinat.sf.sfa import SymmetricFunctionAlgebra_generic
from sage.misc.cachefunc import cached_method
from sage.misc.misc_c import prod
from sage.rings.all import RationalField
from sage.rings.integer_ring import IntegerRing
from sage.rings.semirings.non_negative_integer_semiring import NonNegativeIntegerSemiring

def free_group_elm_to_partition(elm):
    r"""
    Given an element of an abelian free group indexed by natural number
    return an integer sequence of the powers.
    """
    power_dict = elm.dict()
    max_nonzero_entry = max(power_dict.keys())
    return tuple(power_dict.get(i,0) for i in range(max_nonzero_entry+1))

def shifting_operator_action_algebra_elm_to_partition_list(elm):
    elm_list = list(elm)
    return [(free_group_elm_to_partition(supp),coeff) for (supp,coeff) in elm_list]

class ShiftingOperatorActionAlgebra(GroupAlgebra_class):
    r"""
    This is an implementation of the ring on which shifting 
    operators formally act, and is meant to be isomorphic to
    `R[x_1^\pm, x_2^\pm, x_3^\pm, \ldots]` for some ring `R`.
    In particular, the partition `\lambda=(\lambda_1,\lambda_2,\ldots,\lambda_\ell)` 
    is encoded as `x_1^{\lambda_1} x_2^{\lambda_2} \cdots x_\ell^{\lambda_\ell}`
    and this notion generalizes for any sequence of integers with
    finite support.

    Then, the raising operator `R_{ij}` is encoded in this space as simply `\frac{x_i}{x_j}`. 
    This is consisten with how many references work formally with raising operators
    """
    def __init__(self, base_ring, prefix='F'):
        F = Groups.Commutative().free(NonNegativeIntegerSemiring(),prefix)
        category = F.category().Algebras(base_ring)

        GroupAlgebra_class.__init__(self, base_ring, F, category=category)
#        self._register_coercions()
        base = self.base_ring()
        sym = SymmetricFunctions(base)
        h = sym.h()
        h_mor = self.module_morphism(lambda supp: self._supp_to_h(supp,h), codomain=h)
        s = sym.s()
        s_mor = self.module_morphism(lambda supp: self._supp_to_s(supp,s), codomain=s)
        self._outgoing_morphisms = {
            h:h_mor, s:s_mor
        }
    def _element_constructor_(self, seq):
        F = self.group()
        return self.basis()[prod([F.gens()[i]**(seq[i]) for i in range(len(seq))],F.one())]

    def _repr_(self):
        return "Ring of polynomials in countably infinite variables on which shifting operators act over " + repr(self.base_ring())

    def from_iterable_indexed_parent(self, iterable_indexed_elm):
        R = iterable_indexed_elm.parent()
        if self.base_ring().has_coerce_map_from(R.base_ring()):
            f = R.module_morphism(self._element_constructor_, codomain=self)
            return f(iterable_indexed_elm)
        else:
            return self(iterable_indexed_elm)
            
    # def _coerce_map_from_(self, R):
    #     if isinstance(R, SymmetricFunctionAlgebra_generic) or isinstance(R, ShiftingOperatorAlgebra):
    #         if self.base_ring().has_coerce_map_from(R.base_ring()):
    #             return R.module_morphism(self._element_constructor_, codomain=self)
    #         return None
    #     return GroupAlgebra_class._coerce_map_from_(self, R)

    def register_outgoing_morphism(self, support_map, codomain):
        self._outgoing_morphisms[codomain] = self.module_morphism(support_map, codomain=codomain)

    def has_outgoing_morphism(self, codomain):
        return codomain in self._outgoing_morphisms
    
    def outgoing_morphism(self, codomain):
        return self._outgoing_morphisms.get(codomain, None)
        
    def _supp_to_h(self, supp, basis):
        r"""
        Given the support of an element `x_1^{\gamma_1} x_2^{\gamma_2} \cdots x_\ell^{\gamma_\ell}` in the `` ShiftingOperatorActionAlgebra`` and a symmetric function algebra basis `b` generated by `\{b_1, b_2, b_3, \ldots\}`, return the element `b_{\gamma_1} b_{\gamma_2} \cdots b_{\gamma_\ell}` where `b_0 = 1` and `b_{-n} = 0` for all positive integers `n`.
        """
        gamma = list(reversed(sorted(free_group_elm_to_partition(supp))))
        if gamma in Partitions():
            return basis(gamma)
        else:
            return basis.zero()

    def _supp_to_s(self, supp, basis):
        r"""
        Referred to as "Schur straightening" in [BMPS], ...
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
        
        gamma = list(free_group_elm_to_partition(supp))
        rho = list(range(len(gamma) - 1, -1, -1))
        combined = [g + r for g, r in zip(gamma, rho)]
        if has_distinct_parts(combined) and has_nonnegative_parts(combined):
            sign = (-1)**number_of_noninversions(combined)
            sort_combined = reversed(sorted(combined))
            new_gamma = [sc - r for sc, r in zip(sort_combined, rho)]
            return sign * basis(Partition(new_gamma))
        else:
            return basis.zero()

    # def _register_coercions(self):
    #     base = self.base_ring()
    #     sym = SymmetricFunctions(base)
    #     h = sym.h()
    #     h.register_coercion(self.module_morphism(lambda supp: self._supp_to_h(supp,h), codomain=h))
    #     s = sym.s()
    #     s.register_coercion(self.module_morphism(lambda supp: self._supp_to_s(supp,s), codomain=s))
        
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

    - ``base_ring`` -- (default ``RationalField()['t']``) the ring you will use on the raising operators.

    - ``prefix`` -- (default ``"R"``) the label for the raising operators.

    EXAMPLES::

        sage: S = ShiftingOperatorAlgebra()
        sage: s = SymmetricFunctions(RationalField()['t']).s()
        sage: h = SymmetricFunctions(RationalField()['t']).h()

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
    def __init__(self, base_ring=RationalField()['t'], prefix='S', basis_indices=ShiftingSequenceSpace()):
        self._prefix = prefix
        self._base_ring = base_ring
        # a single basis index looks like (1, 0, -1, 2), for example
        self._basis_indices = basis_indices
        # category
        category = Algebras(self._base_ring.category()).WithBasis()
        category = category.or_subcategory(category)
        # init
        self._action_algebra = ShiftingOperatorActionAlgebra(base_ring)
        CombinatorialFreeModule.__init__(
            self,
            self._base_ring,
            self._basis_indices,
            category=category,
            prefix=self._prefix,
            bracket=False)
        self._lift_to_action_algebra = self.module_morphism(lambda supp: self._action_algebra(supp), codomain=self._action_algebra)

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
        elm1 = self._action_algebra(index1)
        elm2 = self._action_algebra(index2)
        res = free_group_elm_to_partition((elm1*elm2).support_of_term()) 
        return self.__getitem__(res)
        # # pad with 0's
        # max_len = max(len(index1), len(index2))
        # index1 = tuple(index1) + (0,) * (max_len - len(index1))
        # index2 = tuple(index2) + (0,) * (max_len - len(index2))
        # # add the vectors
        # index_product = tuple(i1 + i2 for i1, i2 in zip(index1, index2))
        # return self.__getitem__(index_product)

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

        # def _call_basis_on_index(self, seq, index):
        #     """ For internal use only!

        #     Return the action of the basis element indexed by ``seq`` upon the composition ``index``.

        #     INPUTS:

        #     - ``seq`` -- The sequence of the basis element that acts.

        #     - ``index`` -- A sequence (typically a composition or a partition) that we act upon.

        #     EXAMPLES::

        #         sage: S = ShiftingOperatorAlgebra()
        #         sage: S[2, 1]._call_basis_on_index([1, 1], [1, 2, 3, 4, 5])
        #         [2, 3, 3, 4, 5]

        #     ..  SEEALSO::
        #         :meth:`__call__`
        #     """
        #     assert _is_sequence(index)
        #     # pad sequence and index with 0's
        #     index = list(index) + [0] * (len(seq) - len(index))
        #     seq = tuple(seq) + (0,) * (len(index) - len(seq))
        #     # raise and drop
        #     return [v + s for v, s in zip(index, seq)]

        def __call__(self, operand):
            A = self.parent()._action_algebra
            op = A.from_iterable_indexed_parent(self)
            if hasattr(operand,"parent"):
                lift_operand = A.from_iterable_indexed_parent(operand)
                res = op*lift_operand
                operand_parent = operand.parent()
                if A.has_outgoing_morphism(operand_parent):
                    f = A.outgoing_morphism(operand_parent)
                    return f(res)
                else:
                    return shifting_operator_action_algebra_elm_to_partition_list(res)
            else:
                return self(Composition(operand))

                
        # def __call__(self, operand):
        #     r""" Return the action of this shifting operator element on the index ``operand``.

        #     EXAMPLES::

        #         sage: S = ShiftingOperatorAlgebra()

        #         sage: S[2, 1]([1, 2, 3, 4])
        #         [([3, 3, 3, 4], 1)]

        #         sage: S[2, 1](([1, 2, 3, 4], [1, 1, 1]))
        #         ([([3, 3, 3, 4], 1)], [([3, 2, 1], 1)])

        #         sage: Sym = SymmetricFunctions(RationalField()['t'])
        #         sage: s = Sym.schur()
        #         sage: S[2, 1](s[4, 3, 2, 1])
        #         s[6, 4, 2, 1]
        #     """
        #     def raise_func(seq, operand):
        #         # seq is the index
        #         if _is_sequence(operand):
        #             return self._call_basis_on_index(seq, operand)
        #         else:
        #             # it's some symmetric function basis element
        #             parent_basis = operand.parent()
        #             # but does the base ring of the parent basis agree with the base ring of the operator algebra??
        #             # process the vectors
        #             dic = operand.monomial_coefficients()
        #             assert len(dic) == 1
        #             # occasionally a coefficient can show up (not cool, so consider the inclusion of coeff here a patch)
        #             (composition, coeff) = dic.items()[0]
        #             out_composition = raise_func(seq, composition)
        #             if parent_basis.__class__.__name__ in (
        #                     'SymmetricFunctionAlgebra_homogeneous_with_category',
        #                     'SymmetricFunctionAlgebra_elementary_with_category',
        #                     'SymmetricFunctionAlgebra_power_with_category',
        #                     'SymmetricFunctionAlgebra_witt_with_category',
        #                     'SymmetricFunctionAlgebra_schur_with_category',
        #                     'HallLittlewood_qp_with_category'):
        #                 return coeff * straighten(parent_basis, out_composition)
        #             else:
        #                 return coeff * parent_basis(out_composition)

        #     def call_monomial(seq, coeff, operand, power=1):
        #         for _ in range(power):
        #             operand = raise_func(seq, operand)
        #         return (operand, coeff)
        #     # start here
        #     if hasattr(operand, '_get_indices_for_index_operator'):
        #         indices = operand._get_indices_for_index_operator()
        #         # TODO: see if this works for TUPLES of indices.  This has only been tested for a single index.
        #         new_indices = self.__call__(indices)[0][0]
        #         out = operand._new_object_for_index_operator(new_indices)
        #         return out
        #     elif isinstance(operand, tuple):
        #         # the operand is actually a tuple of operands, so perform __call__ on each piece
        #         return tuple(self.__call__(op) for op in operand)
        #     elif _is_sequence(operand):
        #         # the operand is some kind of composition
        #         return [call_monomial(index, coeff, operand) for index, coeff in self]
        #     else:
        #         # the operand is a symmetric function
        #         if len(operand) > 1:
        #             # the operand looks like s[2, 1] + s[3], for example
        #             return sum(self.__call__(summand) for summand in operand.terms())
        #         else:
        #             out_list = [call_monomial(index, coeff, operand)
        #                         for index, coeff in self]
        #             return sum(coeff * mon for mon, coeff in out_list)


class RaisingOperatorAlgebra(ShiftingOperatorAlgebra):
    r""" An algebra of raising operators.

    This class subclasses :class:`ShiftingOperatorAlgebra` and inherits the large majority of its functionality from there.

    We follow the following convention!:

    R[(1, 0, -1)] is the raising operator that raises the first part by 1 and lowers the third part by 1.

    For a definition of raising operators, see [cat]_ Definition 2.1, but be wary that the notation is different there.  See :meth:`ij` for a way to create operators using the notation in the paper.

    If you do NOT want any restrictions on the allowed sequences, use :class:`ShiftingOperatorAlgebra` instead of :class:`RaisingOperatorAlgebra`.

    OPTIONAL ARGUMENTS:

    - ``base_ring`` -- (default ``RationalField()['t']``) the ring you will use on the raising operators.

    - ``prefix`` -- (default ``"R"``) the label for the raising operators.

    EXAMPLES::

        sage: R = RaisingOperatorAlgebra()
        sage: s = SymmetricFunctions(RationalField()['t']).s()
        sage: h = SymmetricFunctions(RationalField()['t']).h()

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

    def __init__(self, base_ring=RationalField()['t'], prefix='R'):
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
        
