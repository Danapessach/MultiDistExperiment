from optparse import OptionParser
import viff.boost
viff.boost.install()
from twisted.internet import reactor

from viff.field import GF
from viff.runtime import Runtime, create_runtime, gather_shares
from viff.config import load_config
from viff.util import rand, find_prime
#from viff import comparison

# We start by defining the protocol, it will be started at the bottom
# of the file.

class Protocol:

    def __init__(self, runtime):
        # Save the Runtime for later use
        self.runtime = runtime

        # This is the value we will use in the protocol.
        print "I am Party %d and my value is %d." \
            % (runtime.id, input)

        # For the comparison protocol to work, we need a field modulus
        # bigger than 2**(l+1) + 2**(l+k+1), where the bit length of
        # the input numbers is l and k is the security parameter.
        # Further more, the prime must be a Blum prime (a prime p such
        # that p % 4 == 3 holds). The find_prime function lets us find
        # a suitable prime.
        l = runtime.options.bit_length
        k = runtime.options.security_parameter
        Zp = GF(find_prime(2**(l + 1) + 2**(l + k + 1), blum=True))

        # We must secret share our input with the other parties. They
        # will do the same and we end up with three variables
		
        m1, m2, m3, m4, m5,m6 = runtime.shamir_share([1, 2, 3, 4, 5,6], Zp, input)
		
        # Now that everybody has secret shared their inputs we can
        # sum and compare them.
		
        sum = m1+ m2+ m3 + m4 + m5+m6
		
        if comp_type=="ge":
            comp=(sum>=comp_value)
        else:
            comp=(sum<=comp_value)
		
        # The results are secret shared, so we must open them before
        # we can do anything usefull with them.
#        open_sum = runtime.open(sum)
		
#        if comp_type=="ge":
#            comp=(open_sum>=comp_value)
#        else:
#            comp=(open_sum<=comp_value)		

        open_comp = runtime.open(comp)		

        # We will now gather the results and call the
        # self.results_ready method when they have all been received.
        #results = gather_shares([open_sum, open_comp])
        results = gather_shares([open_comp])
        results.addCallback(self.results_ready)

        # We can add more callbacks to the callback chain in results.
        # These are called in sequence when self.results_ready is
        # finished. The first callback acts like a barrier and makes
        # all players wait on each other.
        #
        # The callbacks are always called with an argument equal to
        # the return value of the preceeding callback. We do not need
        # the argument (which is None since self.results_ready does
        # not return anything), so we throw it away using a lambda
        # expressions which ignores its first argument.
        runtime.schedule_callback(results, lambda _: runtime.synchronize())
        # The next callback shuts the runtime down, killing the
        # connections between the players.
        runtime.schedule_callback(results, lambda _: runtime.shutdown())

    def results_ready(self, results):
        # Since this method is called as a callback above, the results
        # variable will contain actual field elements, not just
        # Shares. That makes it very easy to work on them.

        # Let us start by unpacking the list.
        #sum = results[0]
        #comp = results[1]
        comp = results[0]

        print "  comp: {%d}" % comp
                

# Parse command line arguments.
parser = OptionParser()
Runtime.add_options(parser)
options, args = parser.parse_args()

id, players = load_config(args[0])
input = int(args[1])
comp_value = int(args[2])
comp_type = args[3] ##ge - grater than equal ##le - less than equal
	
# Create a deferred Runtime and ask it to run our protocol when ready.
pre_runtime = create_runtime(id, players, 1, options)
pre_runtime.addCallback(Protocol)

# Start the Twisted event loop.
reactor.run()