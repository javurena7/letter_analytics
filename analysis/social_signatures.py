import sys
import numpy as np
sys.path.append('..')
import social_signature_comparison as ssc
from analytics import str_to_timestamp
from pandas import DataFrame
from collections import defaultdict
import matplotlib.pyplot as plt; plt.ion()

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
        for j in range(i+1, len(soc_signs)):
            alt_soc = soc_signs[j]
            d_ref = soc_signs[i].get_d_ref(alt_soc.node, alt_soc.valid_sign)
            soc_signs_ref.append(d_ref[0])
    return soc_signs_self, soc_signs_ref

def plot_soc_sign_hist(ss_self, ss_ref, ax, xlabel, ylabel, labels=True):

    ax.hist(ss_self, bins=20, alpha=.6, density=True, label=r'$d_{self}$')
    ax.hist(ss_ref, bins=20, alpha=.6, density=True, label=r'$d_{ref}$')
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    if labels:
        ax.legend(loc=0, fontsize='x-small')
#    ax.set_title(title)


def plot_all_signatures(min_neigh_r, bin_n_r, out_egos, min_bin=3, bin_type='year', savename=''):

    n_vals, valid_nodes = valid_signature_ranges(min_neigh_r, bin_n_r,\
            out_egos, min_bin, bin_type)
    fig, axs = plt.subplots(*n_vals.shape, sharex=True, sharey=True)
    i, j = 0, 0
    for bin_n in bin_n_r:
        print(j)
        for min_neigh in min_neigh_r:
            ss_self, ss_ref = compare_valid_signatures(min_neigh, bin_n, out_egos, valid_nodes, min_bin, bin_type)
            #title = 'min_neigh = {}\n bin_n = {}'.format(min_neigh, bin_n)
            ylabel = 'min neigh=\n{}'.format(min_neigh) if j == 0 else ''
            xlabel = 'bin size={}'.format(bin_n) if (i == n_vals.shape[0] - 1) else ''
            labels = True if (i + j == 0) else False
            plot_soc_sign_hist(ss_self, ss_ref, axs[i,j], xlabel=xlabel, ylabel=ylabel, labels=labels)
            i += 1
        j += 1
        i = 0
    #fig.tight_layout()
    fig.savefig(savename)
    return fig, axs


if __name__ == "__main__":
    min_neigh_r = [5, 7, 10, 15]
    bin_n_r = [1, 3, 5, 10]
    out_egos = out_egos_from_logs(csv_path)
    plot_all_signatures(min_neigh_r, bin_n_r, out_egos, savename='signature_comp_ranges.pdf')

