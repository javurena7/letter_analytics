from ego import *


class CompareSignatureNew(SocialSignature):
    def __init__(self, node, egos, min_neigh, min_bins=3, bin_type='linear', bin_n=3, max_rank=20):
        """
        min_bins: minimium number of consecutive bins necessary to be a valid ego
        min_neigh: minimum number of out-neigbors in a bin to be a valid ego
        """
        super().__init__(node, egos, bin_type, bin_n, max_rank, ns=False)
        self.min_neigh = min_neigh
        self.min_bins = min_bins
        self.valid_sign = None
        self._check_validity()
        self.d_ref = {}


    def _check_validity(self):
        """
        Check if the social signature has at least min_bins consecutive signatures of at least min_neigb valid neighbors.
        """
        cnt = np.array([len(v) for v in self.t_signatures.values()])
        cnt_n = len(cnt) - self.min_bins
        idx = None
        if cnt_n < 0:
            self.valid = False
        else:
            y = [np.all(cnt[i:(i + self.min_bins)] >= self.min_neigh) for i in range(cnt_n)]
            if sum(y) <= 0:
                self.valid = False
            elif sum(y) == 1:
                self.valid = True
                idx = np.where(y)[0][0]
            else:
                # If many sequences are valid, select the one with the most letters
                self.valid = True
                w_list = [self.w_signatures[i] for i in range(len(self.w_signatures))]
                w = [sum(w_list[i:(i + self.min_bins)]) if y[i] else 0 for i in range(cnt_n)]
                idx = np.argmax(w)
        self._valid_idx = idx


    def valid_sign(self):
        if self.valid_sign is not None:
            return self.valid_sign
        elif self.valid:
            return self.get_valid_sign()
        else:
            return None


    def get_valid_sign(self):
        sign = [self.t_signatures[i] for i in range(self._valid_idx, self._valid_idx + self.min_bins)]
        #edges = self.get_year_edges()[self._valid_idx: self._valid_idx + self.min_bins + 1]
        #self.valid_edges = edges
        self.valid_sign = sign
        return sign


    def get_d_self(self):
        """
        Compute self distance between time intervals
        """
        d_t = []
        for s0, s1 in zip(self.valid_sign[:-1], self.valid_sign[1:]):
            p = [s for s in s0.values()]
            q = [s for s in s1.values()]
            d_t.append(self.js(p, q))
        self.d_self = d_t
        self.d_self_avg = np.mean(d_t)


    def add_d_ref(self, node, d_ref):
        done = self.d_ref.get(node, None)


    def get_d_ref(self, node, valid_signs):
        d_t = []
        for s0, s1 in zip(self.valid_sign, valid_signs):
            p = [s for s in s0.values()]
            q = [s for s in s1.values()]
            d_t.append(self.js(p, q))
        self.d_ref[node] = (np.mean(d_t), d_t)
        return self.d_ref[node]


    def update_compare(self, min_neigh, min_bins, bin_type=None, bin_n=None):
        """
        Update number of min neighs and min bins
        New binning only created if both bin_type and bin_n are provided
        """
        if bin_type is not None and bin_n is not None:
            self.update(bin_type, bin_n)
        self.min_neigh = min_neigh
        self.min_bins
        self._check_validity()


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



