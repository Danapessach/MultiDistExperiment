# -*- coding: utf-8 -*-
#
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

"""VIFF runtime. This is where the virtual ideal functionality is
hiding! The runtime is responsible for sharing inputs, handling
communication, and running the calculations.

Each player participating in the protocol will instantiate a
:class:`Runtime` object and use it for the calculations.

The Runtime returns :class:`Share` objects for most operations, and
these can be added, subtracted, and multiplied as normal thanks to
overloaded arithmetic operators. The runtime will take care of
scheduling things correctly behind the scenes.
"""

import time, operator
import struct
from optparse import OptionParser, OptionGroup
from collections import deque
import os
import sys

from viff.field import GF,FieldElement
from viff.util import wrapper, rand
import viff.reactor

from twisted.internet import reactor
from twisted.internet.task import LoopingCall
from twisted.internet.error import ConnectionDone, CannotListenError
from twisted.internet.defer import Deferred, DeferredList, gatherResults
from twisted.internet.defer import maybeDeferred
from twisted.internet.protocol import ReconnectingClientFactory, ServerFactory
from twisted.protocols.basic import Int16StringReceiver


class Share(Deferred):
    """A shared number.

    The :class:`Runtime` operates on shares, represented by this class.
    Shares are asynchronous in the sense that they promise to attain a
    value at some point in the future.

    Shares overload the arithmetic operations so that ``x = a + b``
    will create a new share *x*, which will eventually contain the
    sum of *a* and *b*. Each share is associated with a
    :class:`Runtime` and the arithmetic operations simply call back to
    that runtime.
    """

    def __init__(self, runtime, field, value=None):
        """Initialize a share.

        If an initial value is given, it will be passed to
        :meth:`callback` right away.
        """
        assert field is not None, "Cannot construct share without a field."
        assert callable(field), "The field is not callable, wrong argument?"

        Deferred.__init__(self)
        self.runtime = runtime
        self.field = field
        if value is not None:
            self.callback(value)

    def __neg__(self):
        """Negation."""
        return self.runtime.neg(self)

    def __add__(self, other):
        """Addition."""
        return self.runtime.add(self, other)

    def __radd__(self, other):
        """Addition (reflected argument version)."""
        return self.runtime.add(self, other)

    def __sub__(self, other):
        """Subtraction."""
        return self.runtime.sub(self, other)

    def __rsub__(self, other):
        """Subtraction (reflected argument version)."""
        return self.runtime.sub(other, self)

    def __mul__(self, other):
        """Multiplication."""
        return self.runtime.mul(self, other)

    def __rmul__(self, other):
        """Multiplication (reflected argument version)."""
        return self.runtime.mul(self, other)

    def __invert__(self):
        """Inversion."""
        return self.runtime.invert(self)

    def __pow__(self, exponent):
        """Exponentation to known integer exponents."""
        return self.runtime.pow(self, exponent)

    def __xor__(self, other):
        """Exclusive-or."""
        return self.runtime.xor(self, other)

    def __rxor__(self, other):
        """Exclusive-or (reflected argument version)."""
        return self.runtime.xor(self, other)

    def __lt__(self, other):
        """Strictly less-than comparison."""
        # self < other <=> not (self >= other)
        return 1 - self.runtime.greater_than_equal(self, other)

    def __le__(self, other):
        """Less-than or equal comparison."""
        # self <= other <=> other >= self
        return self.runtime.greater_than_equal(other, self)

    def __gt__(self, other):
        """Strictly greater-than comparison."""
        # self > other <=> not (other >= self)
        return 1 - self.runtime.greater_than_equal(other, self)

    def __ge__(self, other):
        """Greater-than or equal comparison."""
        # self >= other
        return self.runtime.greater_than_equal(self, other)

    def __eq__(self, other):
        """Equality testing."""
        return self.runtime.equal(self, other)

    def __neq__(self, other):
        """Negated equality testing."""
        return 1 - self.runtime.equal(self, other)


    def clone(self):
        """Clone a share.

        Works like :meth:`util.clone_deferred` except that it returns a new
        :class:`Share` instead of a :class:`Deferred`.
        """

        def split_result(result):
            clone.callback(result)
            return result
        clone = Share(self.runtime, self.field)
        self.addCallback(split_result)
        return clone


class ShareList(Share):
    """Create a share that waits on a number of other shares.

    Roughly modelled after the Twisted :class:`DeferredList`
    class. The advantage of this class is that it is a :class:`Share`
    (not just a :class:`Deferred`) and that it can be made to trigger
    when a certain threshold of the shares are ready. This example
    shows how the :meth:`pprint` callback is triggered when *a* and
    *c* are ready:

    >>> from pprint import pprint
    >>> from viff.field import GF256
    >>> a = Share(None, GF256)
    >>> b = Share(None, GF256)
    >>> c = Share(None, GF256)
    >>> shares = ShareList([a, b, c], threshold=2)
    >>> shares.addCallback(pprint)           # doctest: +ELLIPSIS
    <ShareList at 0x...>
    >>> a.callback(10)
    >>> c.callback(20)
    [(True, 10), None, (True, 20)]

    The :meth:`pprint` function is called with a list of pairs. The first
    component of each pair is a boolean indicating if the callback or
    errback method was called on the corresponding :class:`Share`, and
    the second component is the value given to the callback/errback.

    If a threshold less than the full number of shares is used, some
    of the pairs may be missing and :const:`None` is used instead. In
    the example above the *b* share arrived later than *a* and *c*,
    and so the list contains a :const:`None` on its place.
    """
    def __init__(self, shares, threshold=None):
        """Initialize a share list.

        The list of shares must be non-empty and if a threshold is
        given, it must hold that ``0 < threshold <= len(shares)``. The
        default threshold is ``len(shares)``.
        """
        assert len(shares) > 0, "Cannot create empty ShareList"
        assert threshold is None or 0 < threshold <= len(shares), \
            "Threshold out of range"

        Share.__init__(self, shares[0].runtime, shares[0].field)

        self.results = [None] * len(shares)
        if threshold is None:
            self.missing_shares = len(shares)
        else:
            self.missing_shares = threshold

        for index, share in enumerate(shares):
            share.addCallbacks(self._callback_fired, self._callback_fired,
                               callbackArgs=(index, True),
                               errbackArgs=(index, False))

    def _callback_fired(self, result, index, success):
        self.results[index] = (success, result)
        self.missing_shares -= 1
        if not self.called and self.missing_shares == 0:
            self.callback(self.results)
        return result


class ShareLists(Share):
    """Create a share that waits on a number of other shares in a nested list.

    SH: Similar to ShareList. Not used for openings, so no threshold required.
    """

    def __init__(self, shares):
        """Initialize a share list.

        The list of shares must be non-empty and if a threshold is
        given, it must hold that ``0 < threshold <= len(shares)``. The
        default threshold is ``len(shares)``.
        """

        assert len(shares) > 0, "Cannot create empty ShareList"
        
        shareLists = [None]*len(shares)
        for i in xrange(len(shares)):
            if reduce(operator.or_, map(lambda z: isinstance(z, list), shares[i])):
                shareLists[i] = ShareLists(shares[i])
            else:
                shareLists[i] = ShareList(shares[i])
		
        Deferred.__init__(self)

        self.results = [None] * len(shares)
        self.missing_shareLists = len(shares)

        for index, shareList in enumerate(shareLists):
            shareList.addCallbacks(self._callback_fired, self._callback_fired,
                               callbackArgs=(index, True),
                               errbackArgs=(index, False))

    def _callback_fired(self, result, index, success):
        if isinstance(result[0], tuple):
            self.results[index] = [share for (_,share) in result]
        else:
		    self.results[index] = result
        self.missing_shareLists -= 1
        if not self.called and self.missing_shareLists == 0:
            self.callback(self.results)
        return result


def gather(shares):
    """
    SH: Appply Share Lists to each row of shares 
    """	
    if isinstance(shares, list):
        if reduce(operator.or_, map(lambda z: isinstance(z, list), shares)):
            return ShareLists(shares)
        else:
            return gather_shares(shares)
    else:
        return gather_shares(shares)
    
	
def gather_shares(shares):
    """Gather shares.

    Roughly modelled after the Twisted :meth:`gatherResults`
    function. It takes a list of shares and returns a new
    :class:`Share` which will be triggered with a list of values,
    namely the values from the initial shares:

    >>> from pprint import pprint
    >>> from viff.field import GF256
    >>> a = Share(None, GF256)
    >>> b = Share(None, GF256)
    >>> shares = gather_shares([a, b])
    >>> shares.addCallback(pprint)           # doctest: +ELLIPSIS
    <ShareList at 0x...>
    >>> a.callback(10)
    >>> b.callback(20)
    [10, 20]
    """
    
    def filter_results(results):
        return [share for (_, share) in results]
    share_list = ShareList(shares)
    share_list.addCallback(filter_results)
    return share_list

def gather_shares_2d(shares):
    """Gather shares.
    
    By flattening the list of shares
    
    """
    
    flatten_shares = [shares[i][j] for i in range(len(shares)) for j in range(len(shares[i]))]
	
    def filter_results(results):
        flatten_list = [share for (_, share) in results]
        return [[flatten_list[len(shares[i])*i+j] for j in range(len(shares[i]))] for i in range(len(shares))]
    
    share_list = ShareList(flatten_shares)
    share_list.addCallback(filter_results)
    return share_list	


class ShareExchanger(Int16StringReceiver):
    """Send and receive shares.

    All players are connected by pair-wise connections and this
    Twisted protocol is one such connection. It is used to send and
    receive shares from one other player.
    """

    def __init__(self):
        self.peer_id = None
        self.lost_connection = Deferred()
        #: Data expected to be received in the future.
        self.incoming_data = {}
        self.waiting_deferreds = {}
        #: Statistics
        self.sent_packets = 0
        self.sent_bytes = 0

    def connectionMade(self):
        self.sendString(str(self.factory.runtime.id))

    def connectionLost(self, reason):
        reason.trap(ConnectionDone)
        self.lost_connection.callback(self)

    recvd = ""

    def dataReceived(self, recd):
        """Convert int prefixed strings into calls to stringReceived.
        """
        self.recvd = self.recvd + recd
        while len(self.recvd) >= self.prefixLength and not self.paused:
            length ,= struct.unpack(
                self.structFormat,self.recvd[:self.prefixLength])
            if length > self.MAX_LENGTH:
                self.lengthLimitExceeded(length)
                return
            if len(self.recvd) < length+self.prefixLength:
                break
            packet = self.recvd[self.prefixLength:length+self.prefixLength]
            self.recvd = self.recvd[length+self.prefixLength:]
            self.stringReceived(packet)

    def stringReceived(self, string):
        """Called when a share is received.

        The string received is unpacked into the program counter, and
        a data part. The data is passed the appropriate Deferred in
        :class:`self.incoming_data`.
        """
        if self.peer_id is None:
            # TODO: Handle ValueError if the string cannot be decoded.
            self.peer_id = int(string)
            try:
                cert = self.transport.getPeerCertificate()
            except AttributeError:
                cert = None
            if cert:
                # The player ID are stored in the serial number of the
                # certificate -- this makes it easy to check that the
                # player is who he claims to be.
                if cert.get_serial_number() != self.peer_id:
                    print "Peer %s claims to be %d, aborting!" \
                        % (cert.get_subject(), self.peer_id)
                    self.transport.loseConnection()
            self.factory.identify_peer(self)
        else:
            try:
                pc_size, data_size = struct.unpack("!HH", string[:4])
                fmt = "!%dI%ds" % (pc_size, data_size)
                unpacked = struct.unpack(fmt, string[4:])

                program_counter = unpacked[:pc_size]
                data = unpacked[-1]

                key = program_counter
                if key in self.waiting_deferreds:
                    deq = self.waiting_deferreds[key]
                    deferred = deq.popleft()
                    if not deq:
                        del self.waiting_deferreds[key]
                    self.factory.runtime.handle_deferred_data(deferred, data)
                else:
                    deq = self.incoming_data.setdefault(key, deque())
                    deq.append(data)
            except struct.error, e:
                self.factory.runtime.abort(self, e)

    def sendData(self, program_counter, data):
        """Send data to the peer.

        The *program_counter* is a tuple of unsigned integers and 
        *data* is a string.

        The data is encoded as follows::

          +---------+-----------+--------+--------------+
          | pc_size | data_size |   pc   |     data     |
          +---------+-----------+--------+--------------+
            2 bytes   2 bytes     varies      varies

        The program counter takes up ``4 * pc_size`` bytes, the data
        takes up ``data_size`` bytes.
        """
        pc_size = len(program_counter)
        data_size = len(data)
        fmt = "!HH%dI%ds" % (pc_size, data_size)
        t = (pc_size, data_size) + program_counter + (data,)
        packet = struct.pack(fmt, *t)
        self.sendString(packet)
        self.sent_packets += 1
        self.sent_bytes += len(packet)

    def sendShare(self, program_counter, share):
        """Send a share.

        The program counter and the share are converted to bytes and
        sent to the peer.
        """
        self.sendData(program_counter, hex(share.value))

    def loseConnection(self):
        """Disconnect this protocol instance."""
        self.transport.loseConnection()

class SelfShareExchanger(ShareExchanger):

    def __init__(self, id, factory):
        ShareExchanger.__init__(self)
        self.peer_id = id
        self.factory = factory

    def stringReceived(self, program_counter, data):
        """Called when a share is received.

        The string received is unpacked into the program counter, and
        a data part. The data is passed the appropriate Deferred in
        :class:`self.incoming_data`.
        """
        try:
            key = program_counter
            if key in self.waiting_deferreds:
                deq = self.waiting_deferreds[key]
                deferred = deq.popleft()
                if not deq:
                    del self.waiting_deferreds[key]
                self.factory.runtime.handle_deferred_data(deferred, data)
            else:
                deq = self.incoming_data.setdefault(key, deque())
                deq.append(data)
        except struct.error, e:
            self.factory.runtime.abort(self, e)

    def sendData(self, program_counter, data):
        """Send data to the self.id."""
        self.stringReceived(program_counter, data)

    def loseConnection(self):
        """Disconnect this protocol instance."""
        self.lost_connection.callback(self)
        return None


class SelfShareExchangerFactory(ReconnectingClientFactory, ServerFactory):
    """Factory for creating SelfShareExchanger protocols."""

    protocol = SelfShareExchanger
    maxDelay = 3
    factor = 1.234567 # About half of the Twisted default

    def __init__(self, runtime):
        """Initialize the factory."""
        self.runtime = runtime

    def identify_peer(self, protocol):
        raise Exception("Is identify_peer necessary?")

    def clientConnectionLost(self, connector, reason):
        reason.trap(ConnectionDone)

class FakeTransport(object):
    def close(self):
        return True

class ShareExchangerFactory(ReconnectingClientFactory, ServerFactory):
    """Factory for creating ShareExchanger protocols."""

    protocol = ShareExchanger
    maxDelay = 3
    factor = 1.234567 # About half of the Twisted default

    def __init__(self, runtime, players, protocols_ready):
        """Initialize the factory."""
        self.runtime = runtime
        self.players = players
        self.needed_protocols = len(players) - 1
        self.protocols_ready = protocols_ready
#BS inserted following two lines:
        if self.needed_protocols == 0:
            self.protocols_ready.callback(self.runtime)

    def identify_peer(self, protocol):
        self.runtime.add_player(self.players[protocol.peer_id], protocol)
        self.needed_protocols -= 1
        if self.needed_protocols == 0:
            self.protocols_ready.callback(self.runtime)

    def clientConnectionLost(self, connector, reason):
        reason.trap(ConnectionDone)


class Runtime:
    """Basic VIFF runtime with no crypto.

    This runtime contains only the most basic operations needed such
    as the program counter, the list of other players, etc.
    """

    @staticmethod
    def add_options(parser):
        group = OptionGroup(parser, "VIFF Runtime Options")
        parser.add_option_group(group)

        group.add_option("-l", "--bit-length", type="int", metavar="L",
                         help=("Maximum bit length of input numbers for "
                               "comparisons."))
        group.add_option("-k", "--security-parameter", type="int", metavar="K",
                         help=("Security parameter. Comparisons will leak "
                               "information with probability 2**-K."))
        group.add_option("--no-ssl", action="store_false", dest="ssl",
                         help="Disable the use of secure SSL connections.")
        group.add_option("--deferred-debug", action="store_true",
                         help="Enable extra debug output for deferreds.")
        group.add_option("--no-socket-retry", action="store_true",
                         default=False, help="Fail rather than keep retrying "
                         "to connect if port is already in use.")
        group.add_option("--host", metavar="HOST:PORT", action="append",
                         help="Override host and port of players as specified "
                         "in the configuration file. You can use this option "
                         "multiple times on the command line; the first will "
                         "override host and port of player 1, the second that "
                         "of player 2, and so forth.")

        parser.set_defaults(bit_length=32,
                            security_parameter=30,
                            ssl=False,
                            deferred_debug=False)

    def __init__(self, player, threshold, options=None):
        """Initialize runtime.

        Initialized a runtime owned by the given, the threshold, and
        optionally a set of options. The runtime has no network
        connections and knows of no other players -- the
        :func:`create_runtime` function should be used instead to
        create a usable runtime.
        """
#BS-del        assert threshold > 0, "Must use a positive threshold."
        #: ID of this player.
        self.id = player.id
        #: Shamir secret sharing threshold.
        self.threshold = threshold

        if options is None:
            parser = OptionParser()
            self.add_options(parser)
            self.options = parser.get_default_values()
        else:
            self.options = options

        if self.options.deferred_debug:
            from twisted.internet import defer
            defer.setDebugging(True)

        #: Current program counter.
        self.program_counter = [0]

        #: Connections to the other players.
        #:
        #: Mapping from from Player ID to :class:`ShareExchanger`
        #: objects.
        self.protocols = {}

        #: Number of known players.
        #:
        #: Equal to ``len(self.players)``, but storing it here is more
        #: direct.
        self.num_players = 0

        #: Information on players.
        #:
        #: Mapping from Player ID to :class:`Player` objects.
        self.players = {}
        # Add ourselves, but with no protocol since we wont be
        # communicating with ourselves.
        protocol = SelfShareExchanger(self.id, SelfShareExchangerFactory(self))
        protocol.transport = FakeTransport()
        self.add_player(player, protocol)

        #: Queue of deferreds and data.
        self.deferred_queue = deque()
        self.complex_deferred_queue = deque()
        #: Counter for calls of activate_reactor().
        self.activation_counter = 0
        #: Record the recursion depth.
        self.depth_counter = 0
        self.max_depth = 0
        #: Recursion depth limit by experiment, including security margin.
        self.depth_limit = int(sys.getrecursionlimit() / 50)
        #: Use deferred queues only if the ViffReactor is running.
        self.using_viff_reactor = isinstance(reactor, viff.reactor.ViffReactor)

    def add_player(self, player, protocol):
        self.players[player.id] = player
        self.num_players = len(self.players)
        self.protocols[player.id] = protocol

    def shutdown(self):
        """Shutdown the runtime.

        All connections are closed and the runtime cannot be used
        again after this has been called.
        """
        print "Synchronizing shutdown...",

        def close_connections(_):
            print "done."
            print "Closing connections...",
            results = [maybeDeferred(self.port.stopListening)]
            for protocol in self.protocols.itervalues():
                results.append(protocol.lost_connection)
                protocol.loseConnection()
            return DeferredList(results)

        def stop_reactor(_):
            print "done."
            print "Stopping reactor...",
            reactor.stop()
            print "done."

        def close_connection(_):
            print "done."
            print "Closing connection..."
            for protocol in self.protocols.itervalues():
                result = protocol.lost_connection
                protocol.loseConnection()
            return result

        sync = self.synchronize()
        if self.num_players>1:
            sync.addCallback(close_connections)
            sync.addCallback(stop_reactor)
        else:
            sync.addCallback(close_connection)
        return sync

    def abort(self, protocol, exc):
        """Abort the execution due to an exception.

        The *protocol* received bad data which resulted in *exc* being
        raised when unpacking.
        """
        print "*** bad data from Player %d: %s" % (protocol.peer_id, exc)
        print "*** aborting!"
        for p in self.protocols.itervalues():
            p.loseConnection()
        reactor.stop()
        print "*** all protocols disconnected"

    def wait_for(self, *vars):
        """Make the runtime wait for the variables given.

        The runtime is shut down when all variables are calculated.
        """
        dl = DeferredList(vars)
        self.schedule_callback(dl, lambda _: self.shutdown())

    def increment_pc(self):
        """Increment the program counter."""
        self.program_counter[-1] += 1

    def fork_pc(self):
        """Fork the program counter."""
        self.program_counter.append(0)

    def unfork_pc(self):
        """Leave a fork of the program counter."""
        self.program_counter.pop()

    def schedule_callback(self, deferred, func, *args, **kwargs):
        """Schedule a callback on a deferred with the correct program
        counter.

        If a callback depends on the current program counter, then use
        this method to schedule it instead of simply calling
        addCallback directly. Simple callbacks that are independent of
        the program counter can still be added directly to the
        Deferred as usual.

        Any extra arguments are passed to the callback as with
        :meth:`addCallback`.
        """
        self.increment_pc()
        saved_pc = self.program_counter[:]

        @wrapper(func)
        def callback_wrapper(*args, **kwargs):
            """Wrapper for a callback which ensures a correct PC."""
            try:
                current_pc = self.program_counter[:]
                self.program_counter[:] = saved_pc
                self.fork_pc()
                return func(*args, **kwargs)
            finally:
                self.program_counter[:] = current_pc

        return deferred.addCallback(callback_wrapper, *args, **kwargs)

    def schedule_complex_callback(self, deferred, func, *args, **kwargs):
        """Schedule a complex callback, i.e. a callback which blocks a
        long time.

        Consider that the deferred is forked, i.e. if the callback returns
        something to be used afterwards, add further callbacks to the returned
        deferred."""

        if not self.using_viff_reactor:
            return self.schedule_callback(deferred, func, *args, **kwargs)

        if isinstance(deferred, Share):
            fork = Share(deferred.runtime, deferred.field)
        else:
            fork = Deferred()

        def queue_callback(result, runtime, fork):
            runtime.complex_deferred_queue.append((fork, result))

        deferred.addCallback(queue_callback, self, fork)
        return self.schedule_callback(fork, func, *args, **kwargs)

    def synchronize(self):
        """Introduce a synchronization point.

        Returns a :class:`Deferred` which will trigger if and when all
        other players have made their calls to :meth:`synchronize`. By
        adding callbacks to the returned :class:`Deferred`, one can
        divide a protocol execution into disjoint phases.
        """
        self.increment_pc()
        shares = [self._exchange_shares(player, GF(2)(0))
                  for player in self.players]
        result = gather_shares(shares)
        result.addCallback(lambda _: None)
        return result

    def _expect_data(self, peer_id, deferred):
        # Convert self.program_counter to a hashable value in order to
        # use it as a key in self.protocols[peer_id].incoming_data.
        pc = tuple(self.program_counter)
        return self._expect_data_with_pc(pc, peer_id, deferred)

    def _expect_data_with_pc(self, pc, peer_id, deferred):
        key = pc

        if key in self.protocols[peer_id].incoming_data:
            # We have already received some data from the other side.
            deq = self.protocols[peer_id].incoming_data[key]
            data = deq.popleft()
            if not deq:
                del self.protocols[peer_id].incoming_data[key]
            deferred.callback(data)
        else:
            # We have not yet received anything from the other side.
            deq = self.protocols[peer_id].waiting_deferreds.setdefault(key, deque())
            deq.append(deferred)

    def _exchange_shares(self, peer_id, field_element):
        """Exchange shares with another player.

        We send the player our share and record a Deferred which will
        trigger when the share from the other side arrives.
        """
        assert isinstance(field_element, FieldElement)

        if peer_id == self.id:
            return Share(self, field_element.field, field_element)
        else:
            share = self._expect_share(peer_id, field_element.field)
            pc = tuple(self.program_counter)
            self.protocols[peer_id].sendShare(pc, field_element)
            return share

    def _expect_share(self, peer_id, field):
        share = Share(self, field)
        share.addCallback(lambda value: field(long(value, 16)))
        self._expect_data(peer_id, share)
        return share

    def input(self, inputters, field, number=None):
        """Input *number* to the computation.

        The players listed in *inputters* must provide an input
        number, everybody will receive a list with :class:`Share`
        objects, one from each *inputter*. If only a single player is
        listed in *inputters*, then a :class:`Share` is given back
        directly.
        """
        raise NotImplementedError

    def output(self, share, receivers=None):
        """Open *share* to *receivers* (defaults to all players).

        Returns a :class:`Share` to players with IDs in *receivers*
        and :const:`None` to the remaining players.
        """
        raise NotImplementedError

    def add(self, share_a, share_b):
        """Secure addition.
        
        At least one of the arguments must be a :class:`Share`, the
        other can be a :class:`~viff.field.FieldElement` or a
        (possible long) Python integer."""
        raise NotImplementedError

    def mul(self, share_a, share_b):
        """Secure multiplication.
        
        At least one of the arguments must be a :class:`Share`, the
        other can be a :class:`~viff.field.FieldElement` or a
        (possible long) Python integer."""
        raise NotImplementedError

    def handle_deferred_data(self, deferred, data):
        """Put deferred and data into the queue if the ViffReactor is running. 
        Otherwise, just execute the callback."""

        if self.using_viff_reactor:
            self.deferred_queue.append((deferred, data))
        else:
            deferred.callback(data)

    def process_deferred_queue(self):
        """Execute the callbacks of the deferreds in the queue.

        If this function is not called via activate_reactor(), also
        complex callbacks are executed."""

        self.process_queue(self.deferred_queue)

        if self.depth_counter == 0:
            self.process_queue(self.complex_deferred_queue)

    def process_queue(self, queue):
        """Execute the callbacks of the deferreds in *queue*."""

        while queue:
            deferred, data = queue.popleft()
            deferred.callback(data)

    def activate_reactor(self):
        """Activate the reactor to do actual communcation.

        This is where the recursion happens."""

        if not self.using_viff_reactor:
            return

        self.activation_counter += 1

        # setting the number to n makes the reactor called 
        # only every n-th time
        if self.activation_counter >= 2:
            self.depth_counter += 1

            if self.depth_counter > self.max_depth:
                # Record the maximal depth reached.
                self.max_depth = self.depth_counter
                if self.depth_counter >= self.depth_limit:
                    print "Recursion depth limit reached."

            if self.depth_counter < self.depth_limit:
                reactor.doIteration(0)

            self.depth_counter -= 1
            self.activation_counter = 0


def make_runtime_class(runtime_class=None, mixins=None):
    """Creates a new runtime class with *runtime_class* as a base
    class mixing in the *mixins*. By default
    :class:`viff.passive.PassiveRuntime` will be used.
    """
    if runtime_class is None:
        # The import is put here because of circular depencencies
        # between viff.runtime and viff.passive.
        from viff.passive import PassiveRuntime
        runtime_class = PassiveRuntime
    if mixins is None:
        return runtime_class
    else:
        # The order is important: we want the most specific classes to
        # go first so that they can override methods from later
        # classes. We must also include at least one new-style class
        # in bases -- we include it last to avoid overriding __init__
        # from the other base classes.
        bases = tuple(mixins) + (runtime_class, object)
        return type("ExtendedRuntime", bases, {})

def create_runtime(id, players, threshold, options=None, runtime_class=None):
    """Create a :class:`Runtime` and connect to the other players.

    This function should be used in normal programs instead of
    instantiating the Runtime directly. This function makes sure that
    the Runtime is correctly connected to the other players.

    The return value is a Deferred which will trigger when the runtime
    is ready. Add your protocol as a callback on this Deferred using
    code like this::

        def protocol(runtime):
            a, b, c = runtime.shamir_share([1, 2, 3], Zp, input)

            a = runtime.open(a)
            b = runtime.open(b)
            c = runtime.open(c)

            dprint("Opened a: %s", a)
            dprint("Opened b: %s", b)
            dprint("Opened c: %s", c)

            runtime.wait_for(a,b,c)

        pre_runtime = create_runtime(id, players, 1)
        pre_runtime.addCallback(protocol)

    This is the general template which VIFF programs should follow.
    Please see the example applications for more examples.

    """
    
    if runtime_class is None:
        # The import is put here because of circular depencencies
        # between viff.runtime and viff.passive.
        from viff.passive import PassiveRuntime
        runtime_class = PassiveRuntime

    if options and options.host:
        for i in range(len(options.host)):
            players[i + 1].host, port_str = options.host[i].rsplit(":")
            players[i + 1].port = int(port_str)

    # This will yield a Runtime when all protocols are connected.
    result = Deferred()

    # Create a runtime that knows about no other players than itself.
    # It will eventually be returned in result when the factory has
    # determined that all needed protocols are ready.
    runtime = runtime_class(players[id], threshold, options)
    factory = ShareExchangerFactory(runtime, players, result)

    if len(players)==1: 
        return result

    port = players[id].port
    runtime.port = None
    delay = 2
    while runtime.port is None:
        # We keep trying to listen on the port, but with an
        # exponentially increasing delay between each attempt.
        try:
            runtime.port = reactor.listenTCP(port, factory)
        except CannotListenError, e:
            if options and options.no_socket_retry:
                raise
            delay *= 1 + rand.random()
            print "Error listening on port %d: %s" % (port, e.socketError[1])
            print "Will try again in %d seconds" % delay
            time.sleep(delay)
    print "Listening on port %d" % port

    for peer_id, player in players.iteritems():
        if peer_id > id:
            print "Will connect to %s" % player
            reactor.connectTCP(player.host, player.port, factory)

    if runtime.using_viff_reactor:
        # Process the deferred queue after every reactor iteration.
        reactor.setLoopCall(runtime.process_deferred_queue)

    return result

