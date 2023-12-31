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

"""Modeling of Galois (finite) fields. The GF function creates classes
which implements Galois (finite) fields of prime order.

All fields work the same: instantiate an object from a field to get
hold of an element of that field. Elements implement the normal
arithmetic one would expect: addition, multiplication, etc.

Defining a field:

>>> Zp = GF(19)

Defining field elements:

>>> x = Zp(10)
>>> y = Zp(15)
>>> z = Zp(1)

Addition and subtraction (with modulo reduction):

>>> x + y
{6}
>>> x - y
{14}

Bitwise xor for field elements:

>>> z ^ z
{0}
>>> z ^ 0
{1}
>>> 1 ^ z
{0}

Exponentiation:

>>> x**3
{12}

Square roots can be found for elements based on GF fields with a Blum
prime modulus (see :func:`GF` for more information):

>>> x.sqrt()
{3}

Field elements from different fields cannot be mixed, you will get a
type error if you try:

>>> Zq = GF(17)
>>> z = Zq(2)
>>> x + z
Traceback (most recent call last):
    ...
TypeError: unsupported operand type(s) for +: 'GFElement' and 'GFElement'

The reason for the slightly confusing error message is that ``x`` and
``z`` are instances of two *different* classes called ``GFElement``.
"""

from gmpy import mpz
#from math import log


class FieldElement(object):
    """Common base class for elements."""

    def __int__(self):
        """Extract integer value from the field element.

        """
        return self.value

    __long__ = __int__

    def splitBS(self):
        """Splits self into bit array LSB first.

        >>> Zp = GF(29)
        >>> Zp(3).split()
        [{1}, {1}, {0}, {0}, {0}]
        >>> Zp(28).split()
        [{0}, {0}, {1}, {1}, {1}]
        """
        length = int(ceil(log(self.modulus,2)))
        result = [0] * length
        temp = self.value
        for i in range(length):
            result[i] = self.field(temp % 2)
            temp = temp // 2
        return result

#: Cached fields.
#:
#: Calls to GF with identical modulus must return the same class
#: (field), so we cache them here. 
_field_cache = {}

def GF(modulus):
    """Generate a Galois (finite) field with the given modulus.

    The modulus must be a prime:

    >>> Z23 = GF(23) # works
    >>> Z10 = GF(10) # not a prime
    Traceback (most recent call last):
        ...
    ValueError: 10 is not a prime

    Please note, that if you wish to calculate square roots, the
    modulus must be a Blum prime (congruent to 3 mod 4):

    >>> Z17 = GF(17) # 17 % 4 == 1, so 17 is no Blum prime
    >>> x = Z17(10)
    >>> x.sqrt()
    Traceback (most recent call last):
        ...
    AssertionError: Cannot compute square root of {10} with modulus 17
    """
    if modulus in _field_cache:
        return _field_cache[modulus]

    if not mpz(modulus).is_prime():
        raise ValueError("%d is not a prime" % modulus)

    # Define a new class representing the field. This class will be
    # returned at the end of the function.
    class GFElement(FieldElement):

        def __init__(self, value):
            self.value = value % self.modulus

        def __add__(self, other):
            """Addition."""
            if not isinstance(other, (GFElement, int, long)):
                return NotImplemented
            try:
                # We can do a quick test using 'is' here since
                # there will only be one class representing this
                # field.
                assert self.field is other.field, "Fields must be identical"
                return GFElement(self.value + other.value)
            except AttributeError:
                return GFElement(self.value + other)

        __radd__ = __add__

        def __sub__(self, other):
            """Subtraction."""
            if not isinstance(other, (GFElement, int, long)):
                return NotImplemented
            try:
                assert self.field is other.field, "Fields must be identical"
                return GFElement(self.value - other.value)
            except AttributeError:
                return GFElement(self.value - other)

        def __rsub__(self, other):
            """Subtraction (reflected argument version)."""
            return GFElement(other - self.value)

        def __xor__(self, other):
            """Xor for bitvalues."""
            if not isinstance(other, (GFElement, int, long)):
                return NotImplemented
            try:
                assert self.field is other.field, "Fields must be identical"
                return GFElement(self.value ^ other.value)
            except AttributeError:
                return GFElement(self.value ^ other)

        def __rxor__(self, other):
            """Xor for bitvalues (reflected argument version)."""
            return GFElement(other ^ self.value)

        def __mul__(self, other):
            """Multiplication."""
            if not isinstance(other, (GFElement, int, long)):
                return NotImplemented
            try:
                assert self.field is other.field, "Fields must be identical"
                return GFElement(self.value * other.value)
            except AttributeError:
                return GFElement(self.value * other)

        __rmul__ = __mul__

        def __pow__(self, exponent):
            """Exponentiation."""
            return GFElement(pow(self.value, exponent, self.modulus))

        def __neg__(self):
            """Negation."""
            return GFElement(-self.value)

        def __invert__(self):
            """Inversion.

            Note that zero cannot be inverted, trying to do so
            will raise a ZeroDivisionError.
            """
            if self.value == 0:
                raise ZeroDivisionError("Cannot invert zero")

            def extended_gcd(a, b):
                """The extended Euclidean algorithm."""
                x = 0
                lastx = 1
                y = 1
                lasty = 0
                while b != 0:
                    quotient = a // b
                    a, b = b, a % b
                    x, lastx = lastx - quotient*x, x
                    y, lasty = lasty - quotient*y, y
                return (lastx, lasty, a)

            inverse = extended_gcd(self.value, self.modulus)[0]
            return GFElement(inverse)

        def __div__(self, other):
            """Division."""
            try:
                assert self.field is other.field, "Fields must be identical"
                return self * ~other
            except AttributeError:
                return self * ~GFElement(other)

        __truediv__ = __div__
        __floordiv__ = __div__

        def __rdiv__(self, other):
            """Division (reflected argument version)."""
            return GFElement(other) / self

        __rtruediv__ = __rdiv__
        __rfloordiv__ = __rdiv__

        def sqrt(self):
            """Square root.

            No attempt is made the to return the positive square root.

            Computing square roots is only possible when the modulus
            is a Blum prime (congruent to 3 mod 4).
            """
            assert self.modulus % 4 == 3, "Cannot compute square " \
                   "root of %s with modulus %s" % (self, self.modulus)

            # Because we assert that the modulus is a Blum prime
            # (congruent to 3 mod 4), there will be no reminder in the
            # division below.
            root = pow(self.value, (self.modulus+1)//4, self.modulus)
            return GFElement(root)

        def bit(self, index):
            """Extract a bit (index is counted from zero)."""
            return (self.value >> index) & 1

        def signed(self):
            """Return a signed integer representation of the value.

            If x > floor(p/2) then subtract p to obtain negative integer.
            """
            if self.value > ((self.modulus-1)/2):
                return self.value - self.modulus
            else:
                return self.value

        def unsigned(self):
            """Return a unsigned representation of the value"""
            return self.value

        def __repr__(self):
            return "{%d}" % self.value
            #return "GFElement(%d)" % self.value

        def __str__(self):
            """Informal string representation.

            This is simply the value enclosed in curly braces.
            """
            return "{%d}" % self.unsigned()

        def __eq__(self, other):
            """Equality test."""
            try:
                assert self.field is other.field, "Fields must be identical"
                return self.value == other.value
            except AttributeError:
                return self.value == other

        def __ne__(self, other):
            """Inequality test."""
            try:
                assert self.field is other.field, "Fields must be identical"
                return self.value != other.value
            except AttributeError:
                return self.value != other

        def __cmp__(self, other):
            """Comparison."""
            try:
                assert self.field is other.field, "Fields must be identical"
                return cmp(self.value, other.value)
            except AttributeError:
                return cmp(self.value, other)

        def __hash__(self):
            """Hash value."""
            return hash((self.field, self.value))

        def __nonzero__(self):
            """Truth value testing.

            Returns False if this element is zero, True otherwise.
            This allows GF elements to be used directly in Boolean
            formula:

            >>> bool(GF256(0))
            False
            >>> bool(GF256(1))
            True
            >>> x = GF256(1)
            >>> not x
            False
            """
            return self.value != 0

    GFElement.modulus = modulus
    GFElement.field = GFElement

    _field_cache[modulus] = GFElement
    return GFElement
