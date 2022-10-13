# Walter Fontana, 2022
"""
THis module creates a signature object, which holds information about the contact map.
"""
import re
import sys


class KappaSignature:
    """
    This object holds the signature of a system.
    Given a Kappa string specifying a signature, it constructs a representation of the form

        self.signature[agent] = {
                                    site: {
                                            bonds: [ other_site.other_agent, ... ],
                                            states: [ state, state, ... ]
                                          },
                                }

        self.site_types = [a.s, ...]
        self.bond_types[(a1, s1), (a2, s2)] = affinity
             bond types ( (a1, s1), (a2, s2) ) are sorted first by agent type and then by site type
        self.init_agents[a] = init_amount (in nM) or '*' (default)
    """

    def __init__(self, sig):
        """
        'sig' is an expression specifying the Kappa signature. Here we parse it and build the signature object
        """
        # ---------------------------------------------------------------------------------------------------
        # sig structures representing the signature; some redundancy here for convenience

        self.signature_string = sig
        self.signature = {}  # the signature data structure
        self.init_agents = {}  # a dictionary declaring the initial concentration of an agent in nM
        self.site_types = []  # a list of site types
        # A dictionary of bond types and affinities; in the signature they can be indicated
        # as w(eak), m(edium), s(trong), def(ault) or a specific magnitude in nM
        self.bond_types = {}
        # the default state of each agent type
        self.default_agent_state = {}

        # ---------------------------------------------------------------------------------------------------
        # build regex'es

        # change these definitions only if you know what you are doing
        self.symbols = r'[_~][a-zA-Z0-9_~+-]+|[a-zA-Z][a-zA-Z0-9_~+-]*'
        # self.mol = r'(?:@\d+)?'              # integer
        self.mol = r'(?:@[0-9]*[.]?[0-9]+)?'   # floating point number
        self.sID = r'x[0-9]+:'
        self.sep = '.'

        agent_name_re = r'(' + self.symbols + self.mol + r')'
        agent_interface_re = r'\(([^()]*)\)'
        # to dissect agents
        self.agent_re = \
            re.compile(r'^' + agent_name_re + agent_interface_re + r'$')
        # to find all agents (hence groups are non-capturing), before dissecting them
        self.agents_re = \
            re.compile(r'(?:' + self.symbols + self.mol + r')' + r'\([^()]*\)')

        site_name_re = r'(' + self.symbols + r')'
        internal_state_re = r'({.*?})'  # we still will need to parse the state expression
        binding_re = r'(\[.*?\])'  # we still will need to parse the bond expression

        # using optional lookahead, since the internal state is optional and there is no prescribed order.
        # (gobble up the string with .*)
        self.site_re = \
            re.compile(r'^' + site_name_re + r'(?=.*?' + internal_state_re + r')?' + r'(?=.*?' + binding_re + r')?.*')

        # ---------------------------------------------------------------------------------------------------
        # parse the signature

        self.parse_signature(sig)
        self.get_site_and_bond_types()
        # generate the default agent states
        self.default_states()

        # consistency check. invaluable.
        consistent, info = self.consistency()
        if not consistent:
            print('\ninconsistent signature!\n')
            print(info)
            sys.exit()

    def parse_signature(self, sig):
        """
        Parse a Kappa signature expression. The expression can be decorated
        with affinity information and initial agent concentrations (in nM).
        The function fills the signature data structure.
        """
        # capture all agents
        match = self.agents_re.findall(sig)

        for agent in match:
            agent_type, interface = self.parse_sig_agent(agent)
            if agent_type in self.signature:
                sys.exit(f'agent {agent_type} is multiply defined')
            else:
                self.signature[agent_type] = interface

    def parse_sig_agent(self, agent_expression):
        """
        Obtain the agent type and its signature
        """
        match = self.agent_re.match(agent_expression)
        if not match:
            sys.exit('Invalid signature declaration <' + agent_expression + '>')

        agent = match.group(1)

        if '@' in agent:  # we got a concentration decoration
            agent_type, init_amount = agent.split('@')
            self.init_agents[agent_type] = init_amount
        else:
            agent_type = agent
            self.init_agents[agent_type] = '*'  # 'kamix' translates this into a default amount

        # parse the agent interface

        iface = match.group(2)
        iface = iface.replace(',', ' ')
        iface = re.sub(r'\] ', ']*', re.sub(r'} ', '}*', iface))

        if iface == '':
            return agent_type, {}

        interface = {}
        iface_sites = iface.split('*')

        for site in iface_sites:
            site = site.strip()
            match = self.site_re.match(site)
            site_name = match.group(1)
            state_signature = ''
            if match.group(2):  # the modification signature; it may be absent, so we need to check
                state_signature = match.group(2)[1:-1]  # remove parentheses

            binding_signature = ''
            if match.group(3):  # there is an explicit binding signature
                binding_signature = match.group(3)[1:-1]  # remove parentheses
            if site_name in interface:
                sys.exit(f'site {site_name} is multiply defined')
            else:
                interface[site_name] = {'states': state_signature.split(), 'bonds': binding_signature.split()}

        return agent_type, interface

    def default_states(self):
        """
        Generates the default state of each agent according to the signature.
        """
        for agent in self.init_agents:

            if not self.signature[agent]:
                expr = f'{agent}()'
            else:
                expr = f'{agent}('
                for site in self.signature[agent]:
                    expr += f'{site}[.]'
                    if self.signature[agent][site]['states']:  # first state is creation default
                        expr += "{" + f'{self.signature[agent][site]["states"][0]}' + "}"
                    expr += ' '
                s = list(expr)
                s[len(s) - 1] = ')'
                expr = ''.join(s)
            self.default_agent_state[agent] = expr

    def get_site_and_bond_types(self):
        """
        Construct the bond_types dictionary and clean bond lists in 'signature' from decorations.
        Construct the site_types list.
        """
        bond_types = {}
        for agent in self.signature:
            for site in self.signature[agent]:
                self.site_types += [''.join([agent, '.', site])]
                bonds = self.signature[agent][site]['bonds']  # list of bonds
                clean_bonds = []
                for bb in bonds:

                    if '$' in bb:
                        b, affinity = bb.split('$')
                    else:
                        b, affinity = bb, 'def'
                    clean_bonds += [b]

                    site2, agent2 = b.split('.')
                    # this standardizes the bond tuple. 'sorted' sorts the list by sorting the 1st component
                    # and resolving ties by sorting the 2nd component; ascending order.
                    (a1, s1), (a2, s2) = sorted([(agent, site), (agent2, site2)])
                    bnd = (''.join([a1, '.', s1]), ''.join([a2, '.', s2]))
                    if bnd not in bond_types:
                        bond_types[bnd] = affinity
                    else:
                        # a non-def affinity overrides a "def" affinity
                        if bond_types[bnd] != affinity and affinity != 'def':
                            sys.exit(f'inconsistent affinity assignment to bond {bnd}')
                        if bond_types[bnd] == 'def' and affinity != 'def':
                            bond_types[bnd] = affinity

                    self.signature[agent][site]['bonds'] = clean_bonds

        self.bond_types = bond_types

    def consistency(self):
        """
        Performs basic signature sanity checks and tries to provide informative error messages.
        """
        consistent = True
        missing_sites_of_bonds = []  # (A.x, b):  A.x of declared bond b is missing
        wrong_stub = []  # (x.A, B.y):  the stub x.A is missing from B.y
        info = ''
        for b in self.bond_types:
            s1, s2 = b
            agent1, site1 = s1.split('.')
            agent2, site2 = s2.split('.')
            is2 = ''.join([site2, '.', agent2])
            is1 = ''.join([site1, '.', agent1])
            if s1 not in self.site_types:
                missing_sites_of_bonds += [(s1, b)]
            if s2 not in self.site_types:
                missing_sites_of_bonds += [(s2, b)]
            if agent1 not in self.signature:
                info += f'agent {agent1} is not declared\n'
            else:
                if site1 not in self.signature[agent1]:
                    info += f'site {site1} is assumed in bond {b}, but not declared in agent {agent1}\n'
                else:
                    if is2 not in self.signature[agent1][site1]['bonds']:
                        wrong_stub += [(is2, s1)]
            if agent2 not in self.signature:
                info += f'agent {agent2} is not declared\n'
            else:
                if site2 not in self.signature[agent2]:
                    info += f'site {site2} is assumed in bond {b}, but not declared in agent {agent2}\n'
                else:
                    if is1 not in self.signature[agent2][site2]['bonds']:
                        wrong_stub += [(is1, s2)]
        if missing_sites_of_bonds:
            consistent = False
            for (site, bond) in missing_sites_of_bonds:
                info += f'site {site} of declared bond {bond} is missing\n'
        if wrong_stub:
            consistent = False
            for (stub, site) in wrong_stub:
                info += f'stub {stub} is not declared for site {site}\n'

        return consistent, info

    def __str__(self, pp_width=40):
        """
        Pretty print the signature.
        """
        info = f"\n{'SIGNATURE '.ljust(70, '-')}\n\n"
        info += 'signature string\n'
        info += f'{self.signature_string}\n\n'
        for agent in self.signature:
            s = f'agent {agent}'
            info += f'{s:>{pp_width}}\n'
            for site in self.signature[agent]:
                if self.signature[agent][site]["states"]:
                    info += f'{site:>{pp_width}}: states -> {self.signature[agent][site]["states"]}\n'
                if self.signature[agent][site]["bonds"]:
                    info += f'{site:>{pp_width}}:  bonds -> {self.signature[agent][site]["bonds"]}\n'
        info += '\n'
        s = f'{len(self.signature)} agents'
        temp = ''
        for a in self.init_agents:
            temp += f'{a} [{self.init_agents[a]}], '
        info += f'{s:>{pp_width}}: {temp[:-2]}\n'
        s = f'{len(self.site_types)} sites'
        info += f'{s:>{pp_width}}: {self.site_types}\n'
        s = f'{len(self.bond_types)} bond types'
        info += f'{s:>{pp_width}}\n'
        for s1, s2 in self.bond_types:
            b = f"{s1}-{''.join([s2.split('.')[1], '.', s2.split('.')[0]])}"
            info += f'{b:>{pp_width}}: {self.bond_types[(s1, s2)]}\n'
        info += '\n'
        return info


# -----------------------------------------------------------------


if __name__ == '__main__':
    # AP = "A@123(p[a1.P$s a2.P a3.P], l[r.A], r[l.A]), P(a1[p.A], a2[p.A], a3[p.A], d[d.P$w])"
    AP = "A@123(p[a1.P$125.37e-9 a2.P a3.P], l[r.A], r[l.A]), P(a1[p.A], a2[p.A], a3[p.A], d[d.P$w])"
    # e1 = "A(x[a.B y.A]{u p}, y[x.A], z{z1 z2 z3}), B(a[x.A y.B] y[a.B] t{1 2})"
    # e1 = "A(x{u p}[a.B y.A], y[x.A], z{z1 z2 z3}), B(a[x.A y.B] y[a.B] t{1 2})"
    # wrong = "A(x[a.B y.A y.C], y[x.A, v.D]{a,b}, z{z1 z2 z3}, f[m.Z]), B(y[x.A y.B] t{1 2}), C(y[x.A]), D(v[y.A])"
    sig = KappaSignature(AP)
    print(sig)
