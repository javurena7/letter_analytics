import sys
import numpy as np
sys.path.append('..')
import social_signature_comparison as ssc
from analytics import str_to_timestamp
from pandas import DataFrame

csv_path = "/home/javier/Documents/aalto/phd/republic_letters/rrltograph/csv/edges_test.csv"

def out_egos_from_logs(path):
    out_egos = {}
    r = open(path, 'r')
    line = r.readline()

    tags = np.array(line.split('|'))
    node_idx = np.where(tags == 'source')[0][0]
    target_idx = np.where(tags == 'target')[0][0]
    date_idx = np.where(tags == 'start\n')[0][0]
    line = r.readline()
    while line:
        l = line.split('|')
        node = l[node_idx]
        target = l[target_idx]
        date = str_to_timestamp(l[date_idx])

        if node in out_egos:
            if target in out_egos[node]:
                out_egos[node][target].append(date)
            else:
                out_egos[node][target] = [date]
        else:
            out_egos[node] = {target: [date]}
        line = r.readline()
    r.close()
    return out_egos

def get_n_valid_signatures(min_bins, min_neigh, bin_n=10, bin_type='year', out_egos={}, soc_signs=[], ):
    n = 0
    valid_nodes = set()
    if out_egos:
        soc_signs = []
        for node, ego in out_egos.items():
            soc_sign = ssc.CompareSignatureNew(node, ego, min_neigh, min_bins, bin_type, bin_n)
            soc_signs.append(soc_sign)
            if soc_sign.valid:
                n += 1
                valid_nodes.add(node)
    elif soc_signs:
        for soc_sign in soc_signs:
            if soc_sign.bin_type == bin_type and soc_sign.bin_n == bin_n:
                soc_sign.update_compare(min_neigh, min_bins)
            else:
                soc_sign.update_compare(min_neigh, min_bins, bin_type, bin_n)
            if soc_sign.valid:
                n += 1
                valid_nodes.add(soc_sign.node)

    return n, soc_signs, valid_nodes


def valid_signature_ranges(min_neigh_r, bin_n_r, out_egos, min_bin=3, bin_type='year'):
    """
    Compute the number of valid singature for paramter combinations of min_neigh and bin_r
    """
    soc_signs = []

    n_vals = DataFrame(np.zeros((len(min_neigh_r), len(bin_n_r))),
            columns=bin_n_r, index=min_neigh_r)
    valid_nodes_total = {i: defaultdict(dict) for i in bin_n_r}
    i, j = 0, 0
    for bin_n in bin_n_r:
        for min_neigh in min_neigh_r:
            n, soc_signs, valid_set = get_n_valid_signatures(min_bin, min_neigh, bin_n, bin_type, out_egos, soc_signs)
            out_egos = {}
            n_vals.iloc[i][j] = n
            valid_nodes_total[bin_n][min_neigh] = valid_set
            i += 1
        j += 1
        i = 0
    return n_vals, valid_nodes_total


def compare_valid_signatures(min_neigh, bin_n, out_egos, valid_nodes_dict, min_bin=3, bin_type='year'):

    soc_signs = []
    soc_signs_self = []
    valid_nodes = valid_nodes_dict[bin_n][min_neigh]
    for k in valid_nodes:
        s = ssc.CompareSignatureNew(k, out_egos[k], min_neigh, min_bin, bin_type, bin_n)
        _ = s.get_valid_sign()
        s.get_d_self()
        soc_signs.append(s)
        soc_signs_self.append(s.d_self_avg)
    soc_signs_ref = []
    for i in range(len(soc_signs)):
        for j in range(i+1, len(out_egos)):
            alt_soc = soc_signs[j]
            d_ref = soc_signs[i].get_d_ref(alt_soc.node, alt_soc.d_self)
            soc_signs_ref.append(d_ref[0])
    return soc_signs_self, soc_signs_ref







