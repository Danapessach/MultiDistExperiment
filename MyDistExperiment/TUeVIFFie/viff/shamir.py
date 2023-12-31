# Copyright 2007, 2008 VIFF Development Team.
#
# This file is part of VIFF, the Virtual Ideal Functionality Framework.
#
# VIFF is free software: you can redistribute it and/or modify it
# under the terms of the GNU Lesser General Public License (LGPL) as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# VIFF is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
# or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General
# Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with VIFF. If not, see <http://www.gnu.org/licenses/>.

import operator
from viff.util import rand

def share(secret, threshold, num_players):
    """Shamir share secret.

    The *threshold* indicates the maximum number of shares that reveal
    nothing about *secret*. The return value is a list of ``(player
    id, share)`` pairs.

    It holds that sharing and recombination cancels each other:

    >>> from field import GF
    >>> Zp = GF(47)
    >>> secret = Zp(42)
    >>> recombine(share(secret, 7, 15)[:8]) == secret
    True

    The threshold can range from zero (for a dummy-sharing):

    >>> share(Zp(10), 0, 5)
    [({1}, {10}), ({2}, {10}), ({3}, {10}), ({4}, {10}), ({5}, {10})]

    up to but not including *num_players*:

    >>> share(Zp(10), 5, 5)
    Traceback (most recent call last):
      ...
    AssertionError: Threshold out of range
    """
    assert threshold >= 0 and threshold < num_players, "Threshold out of range"

    coef = [secret]
    for j in xrange(threshold):
        # TODO: introduce a random() method in FieldElements so that this 
        # wont have to be a long when we are sharing a GMPIntegerFieldElement.
        coef.append(rand.randint(0, long(secret.modulus)-1))

    shares = []
    for i in xrange(1, num_players+1):
        cur_point = secret.field(i)
        cur_share = coef[threshold]
        for j in xrange(threshold-1, -1, -1):
            cur_share = coef[j] + cur_share * cur_point
        shares.append((cur_point, cur_share))

    return shares

#: The recombination vector used by `recombine` depends only on the
#: recombination point and the player IDs of the shares, and so it can
#: be cached for efficiency.
_recombination_vectors = {}

def recombine(shares, x_recomb=0):
    """Recombines list of ``(xi, yi)`` pairs.

    Shares is a list of *threshold* + 1 ``(player id, share)`` pairs.
    Recombination is done in the optional point *x_recomb*.
    """
    xs = tuple([s[0] for s in shares] + [x_recomb])
    try:
        vector = _recombination_vectors[xs]
    except KeyError:
        vector = []
        for i, x_i in enumerate(xs[:-1]):
            factors = [(x_k - x_recomb) / (x_k - x_i)
                       for k, x_k in enumerate(xs[:-1]) if k != i]
            vector.append(reduce(operator.mul, factors, 1))
        _recombination_vectors[xs] = vector
    sum = 0
    for i in xrange(len(shares)):
        sum += shares[i][1] * vector[i]
    return sum

