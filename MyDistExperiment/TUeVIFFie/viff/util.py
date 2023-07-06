# Copyright 2007, 2008, 2009 VIFF Development Team.
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

"""Miscellaneous utility functions. This module contains various
utility functions used in all parts of the VIFF code. The most
important is the :data:`rand` random generator which is seeded with a
known seed each time. Using this generator for all random numbers
ensures that a protocol run can be reproduced at a later time.
"""

import os
import random
from gmpy import mpz

#: Seed for :data:`rand`.
_seed = os.environ.get('VIFF_SEED')

if _seed is None:
    # If the environment variable is not set, then a random seed is
    # chosen.
    _seed = random.randint(0, 10000)
    print 'Seeding random generator with random seed %d' % _seed
    #: Random number generator used by all VIFF code.
    #:
    #: The generator is by default initialized with a random seed,
    #: unless the environment variable :envvar:`VIFF_SEED` is set to a
    #: value, in which case that value is used instead. If
    #: :envvar:`VIFF_SEED` is defined, but empty, then no seed is used
    #: and a protocol cannot be reproduced exactly.
    rand = random.Random(_seed)
elif _seed == '':
    # If it is set, but set to the empty string (VIFF_SEED=), then no
    # seed is used.
    rand = random.SystemRandom()
else:
    # Otherwise use the seed given, which must be an integer.
    _seed = int(_seed)
    print 'Seeding random generator with seed %d' % _seed
    rand = random.Random(_seed)


def wrapper(func):
    """Decorator used for wrapper functions.

    It is important to use this decorator on any wrapper functions in
    order to ensure that they end up with correct :attr:`__name__` and
    :attr:`__doc__` attributes.
    """
    # Return a decorator which does nothing to the function it is
    # asked to decorate, except update the __name__ and __doc__
    # attributes to match the original wrapped function.
    def decorator(f):
        f.__name__ = func.__name__
        f.__doc__ = func.__doc__
        return f
    return decorator

def find_prime(lower_bound, blum=False):
    """Find a prime above a lower bound.

    If a prime is given as the lower bound, then this prime is
    returned:

    >>> find_prime(37)
    37L

    The bound can be a Python expression as a string. This makes it
    easy for users to specify command line arguments that generates
    primes of a particular bit length:

    >>> find_prime("2**100") # 100 bit prime
    1267650600228229401496703205653L

    Blum primes (a prime p such that p % 4 == 3) can be found as well:

    >>> find_prime(12)
    13L
    >>> find_prime(12, blum=True)
    19L

    If the bound is negative, 2 (the smallest prime) is returned:

    >>> find_prime(-100)
    2L
    """
    lower_bound = eval(str(lower_bound), {}, {})
    if lower_bound < 0:
        prime = mpz(2)
    else:
        prime = mpz(lower_bound - 1).next_prime()

    if blum:
        while prime % 4 != 3:
            prime = prime.next_prime()

    return long(prime)


