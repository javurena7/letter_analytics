from social_signature import *


class CompareSignatures(object):

    def __init__(self, nodes=None, min_deg=15, max_rank=20, result_path='signature/', data_path='1500_1990/30/'):

        self.min_deg = min_deg
        self.max_rank = max_rank
        self.result_path = result_path
        self.data_path = data_path
        if nodes is None:
            nodes = self.get_nodes()
        self.nodes = nodes
        self.average_signature = None

    def get_nodes(self):
        nodes = []
        with open(self.data_path + 'degrees.txt') as r:
            row = r.readline()
            while row:
                row = row.split(' ')
                if int(row[1][:-1]) >= self.min_deg:
                    nodes.append(row[0])
                row = r.readline()
        return nodes

    def _get_ego_networks(self, bin_type='year', bin_n=1):
        sgns = {}
        for node in self.nodes:
            ego = read_ego(node)
            ss = SocialSignature(node, ego, bin_type=bin_type, bin_n=bin_n,
                    max_rank=self.max_rank)
            sgns[node] = ss
        self.signatures = sgns


    def get_year_stats(self, ns=[2, 3, 4, 5, 10]):

        def cnames():
            s = ['node', 'letters']
            for i in [1] + ns:
                s.append('l_' + str(i))
                s.append('w_' + str(i))
            return s

        self._get_ego_networks(bin_type='year', bin_n=1)
        sgns_bins = []
        for node, ss in self.signatures.items():
            w = ss.total_letters
            sgns_bins.append([node, w] + self.compute_longest_period(
                ss.w_signatures.values(), ns))
        sgns_bins = pd.DataFrame(sgns_bins, columns=cnames())
        sgns_bins.to_csv(self.result_path + 'year_stats.csv', index=False)
        self.year_stats = sgns


    def compute_longest_period(self, weights, ns):
        """
        Given a list of social signature weights (pct of letters in a given interval),
        returns the length of the largest consecutive interval, and
        pct of letters sent in that period
        """
        if not isinstance(weights, list):
            weights = list(weights)

        consect = self._build_consect(weights, 1)
        bin_periods = self._get_consect_data(consect, weights, 1)
        for n in ns:
            bin_periods += self._aggregate_weights(weights, n)

        return bin_periods

    def _build_consect(self, weights, n):
        consect = [1]
        for w in weights[1:]:
            if w > 1/(3*len(weights)): # or zero
                if consect[-1] == 0:
                    consect.append(1)
                else:
                    consect[-1] += 1
            else:
                consect.append(0)
        return consect


    def _get_consect_data(self, consect, weights, n):
        arg_mx = np.argmax(consect)
        mx = consect[arg_mx]
        mx_loc = sum(x if x > 0 else 1 for x in consect[:arg_mx])
        try:
            mx_w = sum(weights[mx_loc + 1:mx_loc + mx - 1])
        except:
            mx_w = 0
        return [n * (mx-2), mx_w]


    def _aggregate_weights(self, weights, n):
        """
        Given a list of consecutive periods, aggregate them by 'n' to find potential binnings
        """
        weights_n = [sum(weights[i:i+n]) for i in range(0, len(weights), n)]
        consect = self._build_consect(weights_n, n)
        return self._get_consect_data(consect, weights_n, n)


    def get_average_signature(self, nodes, bin_type, bin_n, remove=None):
        """
        Compute signature given a list of nodes that adhere to certain condition (e.g. length > y)
        """
        avg_sgn = np.zeros(self.max_rank)

        for node in nodes:
            ss = self.signatures[node]
            ss.update(bin_type, bin_n, avg=True)
            avg_sgn += ss.average_signature
        avg_sgn /= len(nodes)
        self.average_signature = avg_sgn

    def reference_distances(self, nodes):
        assert self.average_signature is not None, "Run average signature first"
        refs = []
        selves = []
        for node in nodes:
            ss = self.signatures[node]
            ss.add_reference_distance(self.average_signature)

            refs.append(ss.dist_ref)
            selves.append(ss.dist_self)
        self.dists_ref = refs
        self.dists_self = selves




