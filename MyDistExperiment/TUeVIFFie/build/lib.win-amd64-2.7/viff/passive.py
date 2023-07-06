# -*- coding: utf-8 -*-
#
# Copyright 2008 VIFF Development Team.
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

"""Passively secure VIFF runtime."""

import operator, math

from viff import shamir
from viff.runtime import Runtime, Share, ShareList, gather_shares
from viff.prss import prss, prss_double
from viff.field import GF, FieldElement

from twisted.internet.defer import gatherResults


class PassiveRuntime(Runtime):
    """The VIFF runtime.

    The runtime is used for sharing values (:meth:`shamir_share` or
    :meth:`prss_share`) into :class:`~viff.runtime.Share` object and
    opening such shares (:meth:`open`) again. Calculations on shares
    is normally done through overloaded arithmetic operations, but it
    is also possible to call :meth:`add`, :meth:`mul`, etc. directly
    if one prefers.

    Each player in the protocol uses a :class:`~viff.runtime.Runtime`
    object. To create an instance and connect it correctly with the
    other players, please use the :func:`~viff.runtime.create_runtime`
    function instead of instantiating a :class:`~viff.runtime.Runtime`
    directly. The :func:`~viff.runtime.create_runtime` function will
    take care of setting up network connections and return a
    :class:`Deferred` which triggers with the
    :class:`~viff.runtime.Runtime` object when it is ready.
    """

    def __init__(self, player, threshold, options=None):
        """Initialize runtime."""
        Runtime.__init__(self, player, threshold, options)

    def output(self, share, receivers=None, threshold=None):
        return self.open(share, receivers, threshold)

    def open(self, share, receivers=None, threshold=None):
        """Open a secret sharing.

        The *receivers* are the players that will eventually obtain
        the opened result. The default is to let everybody know the
        result. By default the :attr:`threshold` + 1 shares are
        reconstructed, but *threshold* can be used to override this.

        Communication cost: every player sends one share to each
        receiving player.
        """
        assert isinstance(share, Share)
        # all players receive result by default
        if receivers is None:
            receivers = self.players.keys()
        if threshold is None:
            threshold = self.threshold

        def filter_good_shares(results):
            # Filter results, which is a list of (success, share)
            # pairs.
            return [result[1] for result in results
                    if result is not None and result[0]][:threshold+1]

        def recombine(shares):
            assert len(shares) > threshold
            result = ShareList(shares, threshold+1)
            result.addCallback(filter_good_shares)
            result.addCallback(shamir.recombine)
            return result

        def exchange(share):
            # Send share to all receivers.
            for peer_id in receivers:
                if peer_id != self.id:
                    pc = tuple(self.program_counter)
                    self.protocols[peer_id].sendShare(pc, share)
            # Receive and recombine shares if this player is a receiver.
            if self.id in receivers:
                deferreds = []
                for peer_id in self.players:
                    if peer_id == self.id:
                        d = Share(self, share.field, (share.field(peer_id), share))
                    else:
                        d = self._expect_share(peer_id, share.field)
                        d.addCallback(lambda s, peer_id: (s.field(peer_id), s), peer_id)
                    deferreds.append(d)
                return recombine(deferreds)

        result = share.clone()
        self.schedule_callback(result, exchange)

        # do actual communication
        self.activate_reactor()

        if self.id in receivers:
            return result
    
    def mul_public(self, share_a, share_b):
        assert isinstance(share_a, Share), "share_a should be a share"
        assert isinstance(share_b, Share), "share_b should be a share"
        
        result = gather_shares([share_a, share_b])
        result.addCallback(lambda (a, b): a * b)
        return self.open(result, None, 2*self.threshold) 

    def equal_public(self, share_a, share_b):
        """Public equality test.

        Communication cost: 1 opening (ie. each player sends a share)."""
        field = getattr(share_a, "field", getattr(share_b, "field", None))
        if not isinstance(share_a, Share):
            if not isinstance(share_b, Share):
                return share_a == share_b
            # Comparison with constant: just make a default sharing.
            share_a = Share(self, field, share_a)
        if not isinstance(share_b, Share):
            share_b = Share(self, field, share_b)

        r = self.prss_share_random(field)
        result = gather_shares([share_a, share_b, r])

        def sub_mul_open((a,b,r)):
            # we do not reshare after this multiplication, since
            # we just need to open afterwards anyway
            res = Share(self, field, field(r.value*(a.value - b.value)))
            return self.open(res, None, 2*self.threshold)

        result.addCallback(sub_mul_open)
        result.addCallback(lambda x: x == 0)
        return result

#    def equal_zero_public(self, share_a):
#        """Public zero test.
#
#        Communication cost: 1 opening (ie. each player sends a share)."""
#        field = getattr(share_a, "field", None)
#
#        r = self.prss_share_random(field)
#        result = gather_shares([share_a, r])
#
#        def sub_mul_open((a,r)):
#            # we do not reshare after this multiplication, since
#            # we just need to open afterwards anyway
#            res = Share(self, field, r*a)
#            rr = self.open(res, None, 2*self.threshold)
#            return rr.addCallback(lambda x: x == 0)
#
#        result.addCallback(sub_mul_open)
#        return result

    def equal_zero_public(self, share_a):
        r = self.prss_share_random(share_a.field)
        ra = self.open(r*share_a)
        return ra.addCallback(lambda x: x == 0)

    def neg(self, share_a):
        """Negation of shares.

        Communication cost: none.
        """
        result = share_a.clone()
        result.addCallback(lambda a: -a)
        return result

    def add(self, share_a, share_b):
        """Addition of shares.

        Communication cost: none.
        """
        if not isinstance(share_b, Share):
            # Addition with constant. share_a always is a Share by
            # operator overloading in Share. Clone share_a to avoid
            # changing it.
            result = share_a.clone()
            result.addCallback(lambda a, b: b + a, share_b)
            return result
        
        result = gather_shares([share_a, share_b])
        result.addCallback(lambda (a, b): a + b)
        return result

    def sub(self, share_a, share_b):
        """Subtraction of shares.

        Communication cost: none.
        """
        if not isinstance(share_a, Share):
            share_a = Share(self, share_b.field, share_a)
        if not isinstance(share_b, Share):
            share_b = Share(self, share_a.field, share_b)

        result = gather_shares([share_a, share_b])
        result.addCallback(lambda (a, b): a - b)
        return result

    def lin_comb(self, coefficients, shares):
        """Linear combination of shares.

        Communication cost: none. Saves the construction of unnecessary shares
        compared to using add() and mul()."""

        assert len(coefficients) == len(shares), \
            "Number of coefficients and shares should be equal."

        assert all(map(lambda s: not isinstance(s, Share),coefficients)), \
                "Coefficients should not be shares."

        assert all(map(lambda s: isinstance(s, Share),shares)), \
                "Shares should be shares."

        def computation(shares, coefficients):
            return sum(map(operator.mul, coefficients, shares))

        result = gather_shares(shares)
        result.addCallback(computation, coefficients)
        return result

    def bin_comb(self, shares):
        """Binary combination of shares."""

        def computation(shares):
            sum = 0
            for i in xrange(len(shares)-1, -1, -1):
                sum = 2*sum + shares[i]    
            return sum

        result = gather_shares(shares)
        result.addCallback(computation)
        return result

    def mul(self, share_a, share_b):
        """Multiplication of shares.

        Communication cost: 1 Shamir sharing.
        """
        assert isinstance(share_a, Share), \
            "share_a must be a Share."

        if not isinstance(share_b, Share):
            # Local multiplication. share_a always is a Share by
            # operator overloading in Share. We clone share_a first
            # to avoid changing it.
            result = share_a.clone()
            result.addCallback(lambda a: share_b * a)
            return result

        # At this point both share_a and share_b must be Share
        # objects. So we wait on them, multiply and reshare.

        def share_recombine(number):
            shares = shamir.share(number, self.threshold, self.num_players)

            exchanged_shares = []
            for peer_id, share in shares:
                d = self._exchange_shares(peer_id.value, share)
                d.addCallback(lambda share, peer_id: (peer_id, share), peer_id)
                exchanged_shares.append(d)

            # Recombine the first 2t+1 shares.
            result = gather_shares(exchanged_shares[:2*self.threshold+1])
            result.addCallback(shamir.recombine)
            return result

        result = gather_shares([share_a, share_b])
        result.addCallback(lambda (a, b): a * b)
        self.schedule_callback(result, share_recombine)

        # do actual communication
        self.activate_reactor()

        return result

    def invert(self, share_a):
        r = self.prss_share_random(share_a.field)
        ra = self.open(r*share_a)
        return ra.addCallback(lambda ra, r: ~ra*r, r)

    def gauss(self, A, d, b, c):
        """Gaussian elimination A:= A d - b c

        Communication cost: m * n Shamir sharings.
        """
        
        def share_recombine(number):
            shares = shamir.share(number, self.threshold, self.num_players)

            exchanged_shares = []
            for peer_id, share in shares:
                d = self._exchange_shares(peer_id.value, share)
                d.addCallback(lambda share, peer_id: (peer_id, share), peer_id)
                exchanged_shares.append(d)

            # Recombine the first 2t+1 shares.
            result = gather_shares(exchanged_shares[:2*self.threshold+1])
            result.addCallback(shamir.recombine)
            return result

        def computation(x, m):
            n = len(x)/(m+1) - 1
            mn = m*n
            for i in xrange(m):
                for j in xrange(n):
                    x[n*i+j] = share_recombine(x[n*i+j] * x[mn] - x[mn+1+i] * x[mn+1+m+j])
            del x[mn:]
            return x
        
        def peek(s, k):
            return s[k]
            
        def wrap(s, m, n):
            return [[s.clone().addCallback(peek, n*i+j) for j in xrange(n)] for i in xrange(m)]
            
        result = gather_shares(reduce(operator.add, A) + [d] + b + c)
        self.schedule_callback(result, computation, len(A))

        # do actual communication
        self.activate_reactor()
        
        return wrap(result, len(A), len(A[0]))
        #return result

    def scalar_mul(self, a, x):
        """Scalar multiplication x:= a x

        Communication cost: n Shamir sharings.
        """
        
        def share_recombine(number):
            shares = shamir.share(number, self.threshold, self.num_players)

            exchanged_shares = []
            for peer_id, share in shares:
                d = self._exchange_shares(peer_id.value, share)
                d.addCallback(lambda share, peer_id: (peer_id, share), peer_id)
                exchanged_shares.append(d)

            # Recombine the first 2t+1 shares.
            result = gather_shares(exchanged_shares[:2*self.threshold+1])
            result.addCallback(shamir.recombine)
            return result

        def computation(x):
            a = x.pop()
            for i in xrange(len(x)):
                x[i] = share_recombine(x[i]*a)
            return x    
         
        def peek(s, k):
            return s[k]
            
        def wrap(s, n):
            return [s.clone().addCallback(peek, i) for i in xrange(n)]
            
        result = gather_shares(x+[a])
        self.schedule_callback(result, computation)
        
        # do actual communication
        self.activate_reactor()
        
        return wrap(result, len(x))

    def scalar_prod(self, x, y):
        """Scalar multiplication of vectors x and y

        Communication cost: n Shamir sharings.
        """
        
        def share_recombine(number):
            shares = shamir.share(number, self.threshold, self.num_players)

            exchanged_shares = []
            for peer_id, share in shares:
                d = self._exchange_shares(peer_id.value, share)
                d.addCallback(lambda share, peer_id: (peer_id, share), peer_id)
                exchanged_shares.append(d)

            # Recombine the first 2t+1 shares.
            result = gather_shares(exchanged_shares[:2*self.threshold+1])
            result.addCallback(shamir.recombine)
            return result

        def computation(x):
            n = len(x)/2
            for i in xrange(n):
                x[i] = share_recombine(x[i]*x[n+i])
            del x[n:]
            return x    
        
        def peek(s, k):
            return s[k]
            
        def wrap(s, n):
            return [s.clone().addCallback(peek, i) for i in xrange(n)]
            
        result = gather_shares(x+y)
        self.schedule_callback(result, computation)
        
        # do actual communication
        self.activate_reactor()
        
        return wrap(result, len(x))
  
    def matrix_prod(self, A, B):
        """Computing matrix product of A with transposed B
        using only one round of communication."""
        
        def computation(x, ma, mb):
            n = len(x)/(ma+mb)
            C = []
            for ia in xrange(ma):
                for ib in xrange(mb):
                    s = 0
                    for i in xrange(n):
                        s += x[ia*n+i]*x[(ma+ib)*n+i]
                    C.append(share_recombine(s))
            return C

        def share_recombine(number):
            shares = shamir.share(number, self.threshold, self.num_players)

            exchanged_shares = []
            for peer_id, share in shares:
                d = self._exchange_shares(peer_id.value, share)
                d.addCallback(lambda share, peer_id: (peer_id, share), peer_id)
                exchanged_shares.append(d)

            # Recombine the first 2t+1 shares.
            result = gather_shares(exchanged_shares[:2*self.threshold+1])
            result.addCallback(shamir.recombine)
            return result

        def peek(s, k):
            return s[k]
            
        def wrap(s, ma, mb):
            return [[s.clone().addCallback(peek, mb*ia+ib) for ib in xrange(mb)] for ia in xrange(ma)]

        result = gather_shares(reduce(operator.add, A)+reduce(operator.add, B))
        self.schedule_callback(result, computation, len(A), len(B))

        # do actual communication
        self.activate_reactor()
        
        return wrap(result, len(A), len(B))
                            
    def in_prod(self, shares_a, shares_b):
        """Computing inner product of shares_a with shares_b
        using only one round of communication."""
        
        assert len(shares_a) == len(shares_b), \
            "Number of shares_a and shares_b should be equal."

        def computation(share_list):
            n = len(share_list)/2
            s = 0
            for i in xrange(n):
                s += share_list[i] * share_list[n+i] 
            return s

        def share_recombine(number):
            shares = shamir.share(number, self.threshold, self.num_players)

            exchanged_shares = []
            for peer_id, share in shares:
                d = self._exchange_shares(peer_id.value, share)
                d.addCallback(lambda share, peer_id: (peer_id, share), peer_id)
                exchanged_shares.append(d)

            # Recombine the first 2t+1 shares.
            result = gather_shares(exchanged_shares[:2*self.threshold+1])
            result.addCallback(shamir.recombine)
            return result

        result = gather_shares(shares_a+shares_b)
        result.addCallback(computation)
        self.schedule_callback(result, share_recombine)
        
        # do actual communication
        self.activate_reactor()

        return result
        
    def in_prod2(self, shares_a):
        """Computing inner product of shares_a with itself
        using only one round of communication."""
        
        def computation(share_list):
            n = len(share_list)
            s = 0
            for i in xrange(n):
                s += share_list[i]*share_list[i]
            return s

        def share_recombine(number):
            shares = shamir.share(number, self.threshold, self.num_players)

            exchanged_shares = []
            for peer_id, share in shares:
                d = self._exchange_shares(peer_id.value, share)
                d.addCallback(lambda share, peer_id: (peer_id, share), peer_id)
                exchanged_shares.append(d)

            # Recombine the first 2t+1 shares.
            result = gather_shares(exchanged_shares[:2*self.threshold+1])
            result.addCallback(shamir.recombine)
            return result

        result = gather_shares(shares_a)
        result.addCallback(computation)
        self.schedule_callback(result, share_recombine)
        
        # do actual communication
        self.activate_reactor()

        return result        

    def inner_product(self, list_a, list_b):
        """Multiplication of shares.

        Communication cost: 1 Shamir sharing.
        """
        assert len(list_a) == len(list_b), "lists must be equally long"

        field = getattr(list_a[0], "field", None)

        def share_recombine(number):
            shares = shamir.share(number, self.threshold, self.num_players)

            exchanged_shares = []
            for peer_id, share in shares:
                d = self._exchange_shares(peer_id.value, share)
                d.addCallback(lambda share, peer_id: (peer_id, share), peer_id)
                exchanged_shares.append(d)

            # Recombine the first 2t+1 shares.
            result = gather_shares(exchanged_shares[:2*self.threshold+1])
            result.addCallback(shamir.recombine)
            return result

        def inner((la, lb)):
            return sum(map(operator.mul, la, lb))


        def cont(list_a, list_b):
            def swap(lb, list_a):
                return (list_a, lb)

            result = gather_shares(list_b)
            result.addCallback(swap, list_a)
            return result

        result = gather_shares(list_a)
        result.addCallback(cont, list_b)
        result.addCallback(inner)
        self.schedule_callback(result, share_recombine)

        # do actual communication
        self.activate_reactor()
        return result

    def lsb(self, share):
        """ Least Significant Bit Gate [ST06] """
		
        field = share.field
        l = self.options.bit_length
        k = self.options.security_parameter
        ln = (int)(math.log(self.num_players,2)+1)
		
        assert share.field.modulus > 2**(l+1)+2**(l+k+ln+1), "field is too small"
		
        b = self.prss_share_random(field, True)
        r = self.prss_share_random_max(field, 2**(l+k))
        c = self.open(share+b+2*r)
        self.schedule_callback(c, lambda c,b: c.bit(0)^b, b)
        return c

    def pow(self, share, exponent):
        """Exponentation of a share to an integer by square-and-multiply."""

        assert isinstance(exponent, (int, long)), "Exponent must be an integer"
        assert exponent >= 0, "Exponent must be non-negative"

        if exponent == 0:
            return 1
        elif exponent % 2 == 0:
            tmp = share ** (exponent / 2)
            return tmp * tmp
        else:
            return share * (share ** (exponent-1))

    def xor(self, share_a, share_b):
        field = share_a.field
        if not isinstance(share_b, Share):
            if not isinstance(share_b, FieldElement):
                share_b = field(share_b)
            share_b = Share(self, field, share_b)

        return share_a + share_b - 2 * share_a * share_b
        
        
    def greater_than_equal(self, share_a, share_b):
        """Compute ``share_a >= share_b``.

        Both arguments must be shares from the same field. The result
        is a new 0/1 share from the field.
        """

        field = getattr(share_a, "field", getattr(share_b, "field", None))

        smallField = GF(2147483647)  #Mersenne prime 2^31-1, bit length > 30.9759 old: GF(2703455651027L)
        
        l = self.options.bit_length
        k = self.options.security_parameter # assumption k<=30
        
        # TODO: verify asserts are correct...
        assert field.modulus > 2**(l+2) + 2**(l+k), "Field too small"
        assert smallField.modulus > 3 + 3*l, "smallField too small"
        #BS assert ook smallField.modulus  > k + log(n) .... but PRSS has many values added.
        
        r_bits = [self.prss_share_random(smallField, True) for _ in xrange(l)]
        def convert_bit_share(bit, dst_field):
            """Convert a 0/1 share into *dst_field*."""
            dst_share, src_share = self.prss_share_random_double_max(dst_field, bit.field, 2**(k-1))
            #BS bit.field larger than log n * 2**k
            tmp = self.open(bit + src_share)
            tmp.addCallback(lambda i: dst_field(i.value))
            tmp.field = dst_field
            return tmp - dst_share
        r_bitsField = [convert_bit_share(bit, field) for bit in r_bits]
        r_modl = self.bin_comb(r_bitsField)
        r_divl = self.prss_share_random_max(field, 2**(k-1))
        
        z_rmodl = share_a - share_b + 2**l + r_modl
        c = self.open(z_rmodl + 2**l*r_divl)
        def finish(c, field, smallField, r_bits, z_rmodl):
            s_bit = self.prss_share_random(smallField, True)
            s_sign = 1 - 2 * s_bit
            # mask: uniformly random -- should be non-zero, failure prob. 1/2^k
            mask = self.prss_share_random(smallField, False) 
            #BS mask = mask * mask + 1, assuming -1 is in NQR. 
            
            E = [mask]
            sumXORs = 0
            for i in xrange(l-1, -1, -1):
                c_bit = c.bit(i)
                E.append(s_sign + r_bits[i] - c_bit + 3*sumXORs)
                sumXORs += r_bits[i] ^ c_bit
            E.append(s_sign - 1 + 3*sumXORs)
    
            while len(E) > 1:
                h = []
                while len(E) > 1: h.append(E.pop() * E.pop())
                h.extend(E)
                E = h
      
            E[0] = self.open(E[0])
            E[0].addCallback(lambda bit: field(bit.value != 0))
            UF = E[0] ^ convert_bit_share(s_bit, field)
            z_rmodl = self.open(z_rmodl)
            # return  (z - (c%2**l - r%2**l + UF*2**l)) * 2^(-l)
            return (z_rmodl - (c.value%2**l + UF*2**l)) * ~field(2**l)
        self.schedule_callback(c, finish, field, smallField, r_bits, z_rmodl)
        return c

    def prss_key(self):
        """Create unique key for PRSS.

        This increments the program counter and returns it as a tuple.
        Each straight-line program (typically a callback attached to
        some :class:`Deferred`) is executed in a context with unique
        starting program counter. This ensures that consequetive calls
        to PRSS-related methods will use unique program counters.
        """

        # This is called by every function using PRSS, so do it here.
        # If the assertion is not met, things go wrong, i.e. the PRSS
        # functions generate shares with higher degrees than what
        # open() and mul() expect.
        assert self.threshold >= \
               len(self.players) - len(self.players[self.id].keys.keys()[0]), \
               "PRSS functions have higher threshold than the runtime."

        self.increment_pc()
        return tuple(self.program_counter)

    def prss_share_random(self, field, binary=False):
        """Generate shares of a uniformly random element from the field given.

        If binary is True, a 0/1 element is generated. No player
        learns the value of the element.

        Communication cost: none if binary=False, 1 open otherwise.
        """

        prss_key = self.prss_key()
        prfs = self.players[self.id].prfs(field.modulus)
        share = prss(self.num_players, self.id, field, prfs, prss_key)

        if not binary:
            return Share(self, field, share)

        # Open the square and compute a square-root
        result = self.open(Share(self, field, share*share),
                           threshold=2*self.threshold)
        
        def finish(square, share, binary):
            if square == 0:
                return self.prss_share_random(field, binary)
            else:
                root = square.sqrt()
                # When the root is computed, we divide the share and
                # convert the resulting -1/1 share into a 0/1 share.
                return Share(self, field, (share/root + 1) / 2)

        self.schedule_callback(result, finish, share, binary)
        return result

    def prss_share_random_max(self, field, max):
        prss_key = self.prss_key()
        prfs = self.players[self.id].prfs(max)
        share = prss(self.num_players, self.id, field, prfs, prss_key)
        return Share(self, field, share)

    def prss_share_random_double_max(self, field1, field2, max):
        prss_key = self.prss_key()
        prfs = self.players[self.id].prfs(max)
        share1, share2 = prss_double(self.num_players, self.id, field1, field2, prfs, prss_key)
        return Share(self, field1, share1), Share(self, field2, share2)

    def prss_share_zero(self, field, quantity):
        """Generate *quantity* shares of the zero element from the
        field given.

        Communication cost: none.
        """
        prss_key = self.prss_key()
        prfs = self.players[self.id].prfs(field.modulus)
        zero_share = prss_zero(self.num_players, self.threshold, self.id,
                               field, prfs, prss_key, quantity)
        return [Share(self, field, zero_share[i]) for i in range(quantity)]

    def prss_double_shareBS(self, field, quantity):
        """Make *quantity* double-sharings using PRSS.

        The pair of shares will have degree t and 2t where t is the
        default threshold for the runtime.
        """
        r_t = self.prss_share_random_multi(field, quantity)
        z_2t = self.prss_share_zero(field, quantity)
        return (r_t, [r_t[i] + z_2t[i] for i in range(quantity)])

    def input(self, inputters, field, number=None, threshold=None):
        """Input *number* to the computation.

        The input is shared using the :meth:`shamir_share` method.
        """
        return self.shamir_share(inputters, field, number, threshold)

    def shamir_share(self, inputters, field, number=None, threshold=None):
        """Secret share *number* over *field* using Shamir's method.

        The number is shared using polynomial of degree *threshold*
        (defaults to :attr:`threshold`). Returns a list of shares
        unless there is only one inputter in which case the
        share is returned directly.

        In code it is used like this::

            a, b, c = runtime.shamir_share([1, 2, 3], Zp, x)

        where ``Zp`` is a field and ``x`` is a Python integer holding
        the input of each player (three inputs in total).

        If only a subset of the players provide input it looks like
        this::

            if runtime.id == 1:
                a = runtime.shamir_share([1], Zp, x)
            else:
                a = runtime.shamir_share([1], Zp)

        Instead of branching when calling :meth:`shamir_share`, one
        can give ``None`` as input::

            if runtime.id == 1:
                x = int(raw_input("Input x: "))
            else:
                x = None
            a = runtime.shamir_share([1], Zp, x)

        which might be practical in some cases.

        Communication cost: n elements transmitted.
        """
        assert number is None or self.id in inputters
        if threshold is None:
            threshold = self.threshold

        results = []
        for peer_id in inputters:
            # Unique program counter per input.
            self.increment_pc()

            if peer_id == self.id:
                pc = tuple(self.program_counter)
                shares = shamir.share(field(number), threshold,
                                      self.num_players)
                for other_id, share in shares:
                    if other_id.value == self.id:
                        results.append(Share(self, share.field, share))
                    else:
                        self.protocols[other_id.value].sendShare(pc, share)
            else:
                results.append(self._expect_share(peer_id, field))

        # do actual communication
        self.activate_reactor()

        # Unpack a singleton list.
        if len(results) == 1:
            return results[0]
        else:
            return results
