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

"""Functions for loading and saving player configurations. Each player
participating in a protocol execution must know some information about
the other players, namely their hostname and port number. The player
also needs to know something about itself, namely the keys used for
pseudo-random secret sharing (PRSS).

The :class:`Player` class encapsulates this information. Generating a
player configuration is done using the :func:`generate_configs`
function. The :file:`generate_config_files.py` script uses that
function to generate a player config and save it in a number of
:file:`.ini` files. Such a :file:`.ini` file can be loaded with the
:func:`load_config` function.
"""

from viff.libs.configobj import ConfigObj
from viff.prss import generate_subsets, PRF
from viff.util import rand

class Player:
    """Wrapper for information about a player in the protocol."""

    def __init__(self, id, host, port, keys=None):
        """Initialize a player."""
        self.id = id
        self.host = host
        self.port = port
        self.keys = keys
        self.prfs_cache = {}
        
    def prfs(self, modulus):
        """Retrieve PRSS PRFs.

        The pseudo-random functions are used when this player is part
        of a pseudo-random secret sharing for sharing an element
        random to all players.

        Return a mapping from player subsets to :class:`viff.prss.PRF`
        instances.
        """
        try:
            return self.prfs_cache[modulus]
        except KeyError:
            self.prfs_cache[modulus] = prfs = {}
            for subset, key in self.keys.iteritems():
                prfs[subset] = PRF(key, modulus)
            return prfs

    def __repr__(self):
        """Simple string representation of the player."""
        return "<Player %d: %s:%d>" % (self.id, self.host, self.port)


def load_config(source):
    """Load a player configuration file.

    Configuration files are simple INI-files containing information
    (hostname and port number) about the other players in the
    protocol.

    One of the players own the config file and for this player
    additional information on PRSS keys is available.

    Returns the owner ID and a mapping of player IDs to
    :class:`Player` instances.
    """

    def s_unstr(str):
        """Convert a string to a subset ID."""
        return frozenset(map(int, str.split()))

    def p_unstr(str):
        """Convert a string to a player ID."""
        return int(str[7:])

    if isinstance(source, ConfigObj):
        config = source
    else:
        config = ConfigObj(source, file_error=True)
    players = {}

    for player in config:
        id = p_unstr(player)
        host = config[player]['host']
        port = int(config[player]['port'])

        if 'prss_keys' in config[player]:
            keys = {}
            for subset in config[player]['prss_keys']:
                keys[s_unstr(subset)] = config[player]['prss_keys'][subset]

            players[id] = Player(id, host, port, keys)

            # ID of player for which this config file was made
            owner_id = id
        else:
            players[id] = Player(id, host, port)

    return owner_id, players


def generate_configs(n, t, #BS paillier=ViffPaillier(1024),
                     addresses=None, prefix=None): #BS , skip_prss=False):
    """Generate player configurations.

    Generates *n* configuration objects with a threshold of *t*. The
    *addresses* is an optional list of ``(host, port)`` pairs and
    *prefix* is a filename prefix. One can avoid generating keys for
    PRSS by setting *skip_prss* to True. This is useful when the
    number of players is large.

    The configurations are returned as :class:`ConfigObj` instances
    and can be saved to disk if desired.

    Returns a mapping from player ID to player configuration.
    """
    players = frozenset(range(1, n+1))

    def generate_key():
        # TODO: is a 40 byte hex string as good as a 20 byte binary
        # string when it is used for SHA1 hashing? It ought to be
        # since they contain the same entropy.

        # A SHA1 hash is 160 bit
        return hex(rand.randint(0, 2**160))

    def s_str(subset):
        """Convert a subset to a string."""
        return " ".join(map(str, subset))

    def p_str(player):
        """Convert a player ID to a string."""
        return "Player " + str(player)

    configs = {}
    for p in players:
        config = ConfigObj(indent_type='  ')
        config.filename = "%s-%d.ini" % (prefix, p)
        config.initial_comment = ['VIFF config file for Player %d' % p]
        config.final_comment = ['', 'End of config', '']
        configs[p] = config

    for p in players:
        if addresses is None:
            host, port = 'no-host', 0
        else:
            host, port = addresses[p-1]

        for player, config in configs.iteritems():
            config[p_str(p)] = dict(host=host, port=port)
            if player == p:
                # Prepare the config file for the keys
                config[p_str(p)]['prss_keys'] = {}
                        
    max_unqualified_subsets = generate_subsets(players, n-t)
    for subset in max_unqualified_subsets:
        key = generate_key()
        for player in subset:
            config = configs[player]
            config[p_str(player)]['prss_keys'][s_str(subset)] = key

    return configs
