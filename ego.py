import numpy as np
import pandas as pd
import matplotlib.pyplot as plt; plt.ion()
from scipy.stats import entropy
from collections import OrderedDict, defaultdict, Counter
from analytics import str_to_timestamp, timestamp_to_date


# TODO: remove this from here, read data from a
logs = '../rrltograph/csv/edges_test.csv'



class EgoAnalytics(object):
    """
    Main view for ego networks. Obtains timelines, and interacts with social signature
    """
    def __init__(self, node, out_ego, in_ego):
        self.node = node
        self.out_ego = out_ego
        self.in_ego = in_ego

    def timeline(self, ego):
        ts = sorted([date for node in ego.values() for date in node])
        return ts

    def plot_timelines(self):
        out_ts = self.timeline(self.out_ego)
        in_ts = self.timeline(self.in_ego)
        fig, ax = plt.subplots(figsize=(10, 4))

        self._basic_timeline_plot(out_ts, ax, 'b', 'Sent')
        self._basic_timeline_plot(in_ts, ax, 'r', 'Received')
        ax.legend(loc=0)
        return fig, ax


    def _basic_timeline_plot(self, ts, ax, col, label):
        y_ts = [timestamp_to_date(t).year for t in ts]
        y_c = Counter(y_ts)
        years = sorted(list(y_c.keys()))
        vals = [y_c[y] for y in years]
        ax.plot(years, vals, c=col, label=label)



class SocialSignature(object):
    def __init__(self, node, out_ego, bin_type=None, bin_n=3, max_rank=20):
        """
        node: focus node
        ego_times: dict of outgoing letters
        ego_in_times: dict of incoming letters
        bin_type: 'linear', 'distribution' or 'year'
        bin_n: number of cuts for linear/distribution edge types, or number of years to cut for
        max_rank: maximum number of alters in the social signature
        """
        self.times = out_ego
        self.node = node
        self.max_rank = max_rank
        self.bin_type = bin_type
        self.bin_n = bin_n

        if bin_type is not None:
            self._get_bin_edges(bin_type, bin_n)
            self._get_signatures()
        else:
            bin_n = None


        self.average_signature = None
        self.jsd = None
        self.ns = None

    def update(self, bin_type, bin_n, sim=False, avg=False):
        """
        Update the social signature with new binning
        bin_type = 'linar', 'distribution' or 'year'
        bin_n = if 'linear' or 'distribution', the number of bins;\
                if 'year', the nunber of years in each bin
        sim: if True, compute neighbor similarity (self.ns and self.jsd)
        avg: if True, compute average social signature (self.average_signature)
        """
        self.bin_type = bin_type
        self.bin_n = bin_n
        self._get_bin_edges(bin_type, bin_n)
        self._get_signatures()
        if avg:
            self.average_social_signature()
        if sim:
            self.neighbor_similarity()

    def add_reference_distance(self, ref):
        assert len(ref) == self.max_rank, "Ref. distribution must be of size 'max_rank'"
        self.reference = ref
        self.dist_ref = self.js(ref, self.average_signature)
        if self.jsd is None:
            dist_self = []
            for i, sgn in self.t_signatures.items():
                ref_i = self.average_social_signature(i)
                dist_self.append(self.js(ref_i, sgn.values()))
            self.dist_self = np.mean(dist_self)
        else:
            self.dist_self = np.mean(self.jsd.diagonal())


    def _create_signature(self, times):
        x = dict((k, len(time)) for k, time in times.items())
        l = sum(x.values())
        signature = OrderedDict((k, v/l) for k, v in sorted(x.items(), key=lambda item: item[1], reverse=True))
        return signature


    def _get_bin_edges(self, bin_type, bin_n):
        """
        Create bins, where social sigunatures will be calculated
        """
        s_times = sorted([t for time in self.times.values() for t in time])
        if bin_type == 'linear':
            self.edges = np.linspace(s_times[0] - 1, s_times[-1], bin_n + 1)
        elif bin_type == 'distribution':
            quantiles = np.linspace(0, 100, bin_n + 1)
            self.edges = np.percentile(s_times, quantiles)
            self.edges[0] =- 1 # Edge correction
        elif bin_type == 'year':
            bin_size = 365.25 *  bin_n
            self.edges = np.arange(int(s_times[0]-1), int(s_times[-1]) + bin_size, bin_size)

    def _get_signatures(self):
        """
        Calculate social signatures for each bin
        """
        t_signatures = {i: {} for i in range(len(self.edges) - 1)}
        w_signatures = {i: 0 for i in range(len(self.edges) - 1)}
        for neigh, ts in self.times.items():
            sgn_ts = np.digitize(ts, self.edges, True) - 1
            for t, sgn in zip(ts, sgn_ts):
                t_signatures[sgn].setdefault(neigh, []).append(t)
                w_signatures[sgn] += 1

        self.total_letters = sum(w_signatures.values())
        w_signatures = {k: w/self.total_letters for k, w in w_signatures.items()}
        t_signatures = {k: self._create_signature(v) for k, v in t_signatures.items()}

        self.w_signatures = w_signatures
        self.t_signatures = t_signatures

    def _entropy(self, p):
        return - sum([pi * np.log(pi) for pi in p])

    def js(self, p, q):
        """
        Jensen-Shannon Divergence for measuring differences between distributions
        """
        if (len(p) == 0) or (len(q)) == 0:
            return np.nan
        if len(p) < self.max_rank:
            p = list(p) + [0.] * (self.max_rank - len(p))
        else:
            p = list(p)[:self.max_rank]
        if len(q) < self.max_rank:
            q = list(q) + [0.] * (self.max_rank - len(q))
        else:
            q = list(q)[:self.max_rank]
        p, q = np.array(p), np.array(q)
        p /= p.sum()
        q /= q.sum()
        m = (p + q) / 2
        jsd = entropy((p + q)/2) - (entropy(p) + entropy(q)) / 2
        return np.sqrt(jsd)

    def neighbor_similarity(self):
        """
        Jaccard Similarity for the set of neighbors at different times. Allowing to \
                check whether the set of neighbors changes at different times.
        """
        n_sgn = len(self.t_signatures)
        ns = np.zeros((n_sgn, n_sgn))
        jsd = ns.copy()
        for i in range(n_sgn):
            ref_i = self.average_social_signature(i)
            js_ii = self.js(ref_i, self.t_signatures[i].values())
            jsd[i, i] = js_ii
            for j in range(i + 1, n_sgn):
                s_ij = self.jaccard_similarity(self.t_signatures[i], self.t_signatures[j])
                js_ij = self.js(self.t_signatures[i].values(), self.t_signatures[j].values())
                ns[i, j] = s_ij
                jsd[i, j] = js_ij

        self.ns = ns
        self.jsd = jsd

    def jaccard_similarity(self, sa, sb):
        if (len(sa) == 0) or (len(sb) == 0):
            return np.nan
        sa = set(list(sa.keys())[:self.max_rank])
        sb = set(list(sb.keys())[:self.max_rank])
        try:
            sab = len(sa.intersection(sb)) / len(sa.union(sb))
        except ZeroDivisionError:
            sab = 0
        return sab

    def average_social_signature(self, remove=None):
        """
        Compute average social signature. If remove is an int, compute the
        avg soc sign without the signature of time 'remove'
        """
        n_sign = len(self.t_signatures)
        av_ss = np.zeros(self.max_rank)
        for k, v in self.t_signatures.items():
            if k != remove:
                lv = len(v)
                if lv >= self.max_rank:
                    hpr = list(v.values())[:self.max_rank]
                if lv < self.max_rank:
                    hpr = list(v.values()) + [0.] * (self.max_rank - lv)
                av_ss += self.w_signatures[k] * np.array(hpr)

        if remove is None:
            self.average_signature = av_ss
        else:
            return (1/self.w_signatures[remove]) * av_ss if \
                    self.w_signatures[remove] > 0 else av_ss


    def plot_signature(self, signature, ax=None):
        if ax is None:
            fig, ax = plt.subplots(1)
        y = list(signature.values())[:self.max_rank]
        ax.plot(y, '.')
        if self.average_signature is not None:
            ax.plot(self.average_signature, alpha=.2, c='grey')
        return ax

    def plot_all_signatures(self):
        fig, ax = plt.subplots(len(self.t_signatures))


def read_ego(node):
    """
    Read ego network of sent letters for node
    """
    ego = defaultdict(list)
    ego_in = defaultdict(list)
    with open(logs, 'r') as r:
        line = r.readline()
        t = line.split('|')
        t0_idx, t1_idx, time_idx = t.index('source'), t.index('target'), t.index('start\n')
        while line:
            t = line.split('|')
            if t[t0_idx] == node:
                ego[t[t1_idx]] += [str_to_timestamp(t[time_idx])]
            if t[t1_idx] == node:
                ego_in[t[t0_idx]] += [str_to_timestamp(t[time_idx])]
            line = r.readline()
        ego = {k: sorted(v) for k, v in ego.items()}
        ego_in = {k: sorted(v) for k, v in ego_in.items()}
    return ego, ego_in


if __name__ == '__main__':
    import argparse
    import json
    from os import path

    parser = argparse.ArgumentParser()
    parser.add_argument('--node', default='cccc8972-7adb-463d-b039-bcee2898b222', type=str) #John Locke
    parser.add_argument('--bin_type', default='year', type=str)
    parser.add_argument('--bin_n', default=2, type=int)
    parser.add_argument('--savepath', default='', type=str)
    pargs = parser.parse_args()
    node, bin_type, bin_n, savepath = pargs.node, pargs.bin_type, pargs.bin_n, pargs.savepath

    ego = read_ego(node)
    a = SocialSignature(node, ego, bin_type=bin_type, bin_n=bin_n)
    a.average_social_signature()
    a.neighbor_similarity()

    if savepath:
        spath = path.join(savepath, 'signatures.json')
        with open(spath, 'w') as f:
            json.dump(a.t_signatures, f)


