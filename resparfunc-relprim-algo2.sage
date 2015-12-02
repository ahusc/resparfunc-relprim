"""
resparfunc-relprim-algo2

A restricted partition function counts the number of possible partitions
of a nonnegative integer, subject to various constraints.
This module is about restricted partition functions of a special kind,
also known as Sylvester's denumerant, for the special case of pairwise
relatively prime numbers.
"""
#***********************************************************************
#       Copyright (C) 2015 Anne Huschitt <Anne.Huschitt@gmail.com>,
#                          Alfred Seiler
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#***********************************************************************

class RestrictedPartitionFunctionRelprim(object):
    r"""
    RestrictedPartitionFunctionRelprim allows to get the number
    of partitions with parts in a given list of
    pairwise relatively prime positive integers
    for arbitrary nonnegative integers very efficiently.

    Once constructed with
    ``RestrictedPartitionFunctionRelprim(parts_in=A)``
    for a specific restriction list ``A``, a
    ``RestrictedPartitionFunctionRelprim`` object can be used to
    calculate the number of partitions for an integer ``t`` using the
    object's method ``number_of_partitions(t)``.
    Note that this gives the same result as
    ``Partitions(t, parts_in=A).cardinality()``, but the once
    constructed object can be reused for different values of ``t`` and
    the values of ``t`` can be considerably larger.

    The representation of a ``RestrictedPartitionFunctionRelprim`` object
    is based on the fact that the restricted partition function `p^{(A)}` for
    `A=(a_{1},a_{2},...,a_{m})`, where the positive integers `a_{k}`,
    `k=1,...,m`, are pairwise relatively prime, can be written as
    `p^{(A)}=q + \sum_{k=1}^{m}f_{a_{k}}` where `q` is a 
    polynomial with rational coefficients and degree `m-1` and 
    `f_{a_{k}}` are `m` periodical functions with period `a_{k}`
    each, which satisfy the condition
    `\sum_{i=0}^{a_{k}-1}f_{a_{k}}(i) = 0`.
    We represent a periodical function `f_{a_{k}}` by an object of type
    ``PeriodObject`` containing its period `a_{k}` and a list of
    function values `f_{a_{k}}(0),...,f_{a_{k}}(a_{k}-1)`.

    The time needed to construct a ``RestrictedPartitionFunctionRelprim``
    object for a specific restriction list ``A`` in seconds is printed when
    the calculation is finished, but this can be suppressed by setting the
    class variable ``RestrictedPartitionFunction.print_time=False``.
    
    EXAMPLES:

    Construct partition function for a specific restriction list
    and use it to get the number of partitions for several numbers::

        sage: RestrictedPartitionFunctionRelprim.print_time=False  # suppress timing for doctest
        sage: pa=RestrictedPartitionFunctionRelprim([4,3,5])
        sage: pa.number_of_partitions(8)
        2
        sage: pa.number_of_partitions(10000)
        834334
        sage: pa.number_of_partitions(10**30)
        8333333333333333333333333333433333333333333333333333333334
        sage: pa.number_of_partitions(10**10**6).ndigits()
        1999998

    Showing the polynomial and the periodical functions which represent
    the partition function for a specific restriction list::

        sage: pa=RestrictedPartitionFunctionRelprim([4,3,5])
        sage: pa.polynomial
        1/120*t^2 + 1/10*t + 191/720
        sage: pa.periodic_functions
        [P=4, values=[5/16, -1/16, -3/16, -1/16],
         P=3, values=[2/9, -1/9, -1/9],
         P=5, values=[1/5, -1/5, -1/5, 1/5, 0]]

    Comparing ``RestrictedPartitionFunctionRelprim().number_of_partitions()``
    with ``Partitions().cardinality()``. The result must be the same::

        sage: A=[4, 5, 21]
        sage: pa=RestrictedPartitionFunctionRelprim(parts_in=A)
        sage: t=500
        sage: pa.number_of_partitions(t) == Partitions(t, parts_in=A).cardinality()
        True
        sage: len([t for t in [100..200] if pa.number_of_partitions(t) != Partitions(t, parts_in=A).cardinality()])
        0

    """

    # class variable to influence time output
    print_time = True

    def __init__(self, parts_in=None):
        """
        Initialize ``self``.

        INPUT:

        - ``parts_in`` -- a list of positive integers

        See documentation of RestrictedPartitionFunctionRelprim for details.
        """
        if not isinstance(parts_in, (list, tuple)):
            raise ValueError('%s is not a list or tuple'%parts_in)
        for act_a in parts_in:
            if not act_a in ZZ or not act_a > 0:
                raise ValueError('%s is not a positive integer'%act_a)
        gcdlist = [[gcd(a, b), i, j] for i, a in enumerate(parts_in)
                    for j, b in enumerate(parts_in)
                    if j > i
                  ]
        for el in gcdlist:
            if el[0] > 1:
                raise ValueError(
                    '{0} and {1} are not coprime (gcd={2})'.format(
                    parts_in[el[1]], parts_in[el[2]], el[0])) 

        start_time = cputime()
        self._parts_in = parts_in
        # periodic functions
        self.periodic_functions = []
        act_parts_list = list(parts_in)
        a = parts_in[0]
        m = len(parts_in)
        for k in range(m):
            fa = PeriodObject(period = a,
                    values = [(a - 1) / a] + [-1/a for _ in range(1, a)])
            for l in range(1, m):
                tu = val = ma = 0
                alp1 = act_parts_list[l]
                for _ in range(a):
                    val += fa[tu]
                    fa[tu] = val
                    ma += val
                    tu = (tu + alp1) % a
                ma = ma / a
                for i in range(a):
                    fa[i] -= ma
            self.periodic_functions.append(fa)
            if k < m - 1:
                # put next a at start of list
                a = act_parts_list[k + 1]
                act_parts_list[k + 1] = act_parts_list[0]
                act_parts_list[0] = a
        # calculate values
        vals = []
        funcvals = number_partitions_brute_force(parts_in);
        for t in range(m):
            val = funcvals[t]
            for k in range(m):
                val -= self.periodic_functions[k][t%parts_in[k]]
            vals.append((t, val))
        # interpolate polynomial
        R = PolynomialRing(QQ, 't')
        self.polynomial = R.lagrange_polynomial(vals)    
        cls = self.__class__
        if cls.print_time:
            print "time: " + str(cputime(start_time))

    def number_of_partitions(self, t):
        r"""
        Get the number of partitions of a nonnegative integer ``t`` with
        parts in a specific restriction list as contained in the current
        object of type ``RestrictedPartitionFunctionRelprim``.

        Note that this gives the same result as
        ``Partitions(t, parts_in=A).cardinality()``, but the once
        constructed object of type ``RestrictedPartitionFunctionRelprim``
        can be used to call ``number_of_partitions`` consecutively for
        different values of ``t`` and the values of ``t`` can be
        considerably larger.

        INPUT:

        - ``t`` -- a nonnegative integer


        OUTPUT:

        - the number of partitions of ``t`` with parts in the restriction
          list as contained in the current object

        """
        if not t in ZZ or not t > 0:
            raise ValueError('%s is not a nonnnegative integer'%t)
        sval = self.polynomial(t)
        for pf in self.periodic_functions:
            sval += pf[t % pf.period]
        assert sval in ZZ
        return ZZ(sval)

    def __repr__(self):
        return (
            "RestrictedPartitionFunctionRelprim for {0!r}".
            format(self._parts_in)
        )


class PeriodObject(object):
    """Holds a period object."""
    def __init__(self, values=None, period=None):
        self.period = period
        if isinstance(values, list) and len(values) == period:
            self._vals = values
        else:
            raise AssertionError('conflicting period={0} and values={1}'.format(period, values))
    def __repr__(self):
        return "P={0.period}, values={0._vals}".format(self)
    def __iter__(self):
        return iter(self._vals)
    def __getitem__(self, pos):
        return self._vals[pos]
    def __setitem__(self, pos, val):
        self._vals[pos] = val


def number_partitions_brute_force(parts_in):
    r"""
    Get the number of partitions with parts in a specific restriction
    list for each number starting from 0 up to (but not including) the
    length of the restriction list, using a brute force method.


    INPUT:

    - ``parts_in`` -- a list of positive integers


    OUTPUT:

    - a list with the same length as ``parts_in`` containing the number
      of partitions for consecutive numbers starting with 0
    """
    limit = len(parts_in)
    funcvals = [1] + [0 for _ in range(limit - 1)]
    for act_a in parts_in:
        if act_a >= limit:
            continue
        new_funcvals = [0 for _ in range(limit)]
        a_part = 0
        while a_part <= limit:
            for t in range(limit):
                rest = t - a_part
                if rest >= 0:
                    new_funcvals[t] += funcvals[rest]
            a_part += act_a
        funcvals = new_funcvals
    return funcvals

