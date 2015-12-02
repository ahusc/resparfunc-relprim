"""
fourier-dedekind-sum

Fourier-Dedekind sums are number-theoretic objects. They are closely
related to the periodical functions which are part of the representation
of the restricted partition function for the special case of pairwise
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

def fourier_dedekind_sum(alist, b):
    r"""
    fourier_dedekind_sum allows to get the Fourier-Dedekind sum
    `s_{n}(a_{1},a_{2},...,a_{d};b)` for a list of positive integers
    `a_{1},a_{2},...,a_{d}` and a positive integer `b`, where `b` is
    relatively prime with each `a_{i}`, `i=1,...,m`.

    INPUT:

    - ``alist`` -- a list of positive integers
    - ``b`` -- a positive integer coprime with each integer from `alist`

    OUTPUT:

    - a list of `b` rational numbers representing the function values of
      `s_{n}(a_{1},a_{2},...,a_{d};b)` on a period `0,...,b-1`

    EXAMPLES:
    
    Calculate a small Fourier-Dedekind sum::

        sage: fourier_dedekind_sum([3,5],7)
        [1/7, 1/7, 0, -2/7, 2/7, -2/7, 0]

    Obtain the Fourier-Dedekind sum for a list consisting of the integer
    5 repeated 10 times to illustrate that the algorithm is also usable
    for numbers with are not coprime. Only `b` must be coprime with each
    number of the list. As the resulting Fourier-Dedkind sum is large,
    we do not display it directly, but only the first, the maximum and
    the minimum value::

        sage: fds=fourier_dedekind_sum([5,5,5,5,5,5,5,5,5,5],101)
        sage: fds[0]
        -2198653249650/101
        sage: max(fds)
        2304180724030/101
        sage: min(fds)
        -2309826017050/101

    """

    # parameter validation
    if not b in ZZ or not b > 0:
        raise ValueError('%s is not a positive integer'%b)
    if not isinstance(alist, (list, tuple)):
        raise ValueError('%s is not a list or tuple'%alist)
    for act_a in alist:
        if not act_a in ZZ or not act_a > 0:
            raise ValueError('%s is not a positive integer'%act_a)
        act_gcd = gcd(act_a, b)
        if act_gcd > 1:
            raise ValueError(
                '{0} and {1} are not coprime (gcd={2})'.format(
                act_a, b, act_gcd)) 
    # initialization
    fds = [(b - 1)/b] + [-1/b for _ in range(b-1)]
    # consider each number from list
    for a in alist:
        pos = 0
        newval = fds[0]
        cumsum = 0
        # calculate the b function values in specified sequence
        for _ in range (b):
            pos = (pos - a) % b
            newval = newval + fds[pos]
            fds[pos] = newval
            cumsum = cumsum + newval
        # subtract mean value from each function value
        meanval = cumsum / b
        for s in range(b):
            fds[s] = fds[s] - meanval
    return fds

