import os
import utils
import numpy as np
import pandas as pd
from netpython import *
import datetime as dt
from verkko.binner.binhelp import *
import matplotlib.pyplot as plt
from scipy.stats import rankdata
# Functions for analysis

########################
# percolation
#######################

def lineToEdge(line):
    s,d,w = line.split()
    s,d,w = int(s), int(d), float(w)
    return (s, d, w)


def write_sorted_edges(net, output_path, reverse=True, func=None):
    """
    Sort edges of a net and write a file with their sorted (or not) edges.
    If reverse=True, the biggest weights will go first.
    """
    edges = [edge for edge in net.edges]
    if func:
        edges.sort(key=lambda x: func(x[2]), reverse=reverse)
    else:
        edges.sort(key=lambda x: x[2], reverse=reverse)
    f = open(output_path, "w+")
    for edge in edges:
        e = str(edge[0]) + " " + str(edge[1]) + " " + str(edge[2]) + "\n"
        f.write(e)
    f.close()

def write_sorted_edges_from_dic(dic, output_path, reverse=True):

    edges = [(key[0], key[1], value) for key, value in dic.iteritems()]
    n_edges = len(edges)
    ranks = rankdata([e[2] for e in edges], 'max')/float(n_edges)
    edge_rank = [(edges[i][0], edges[i][1], ranks[i]) for i in range(n_edges)]
    edge_rank.sort(key=lambda x: x[2], reverse=reverse)
    f = open(output_path, "w+")
    for edge in edge_rank:
        e = str(edge[0]) + " " + str(edge[1]) + " " + str(edge[2]) + "\n"
        f.write(e)
    f.close()

def net_overlap(net, output_path=None, alt_net_path=None):
    #if not return_net:
    #    dic = {}
    #    for edge in net.edges:
    #        e0, e1, w = edge
    #        e0, e1 = np.sort([e0, e1])
    #        dic[(e0, e1)] = netanalysis.overlap(net, e0, e1)
    #    return dic
    if output_path is not None:
        f = open(output_path, 'w+')
        if alt_net_path is not None:
            r = open(alt_net_path, 'r')
            row = r.readline()
            while row:
                e0, e1, _ = utils.parse_time_line(row)
                e0, e1 = np.sort([int(e0), int(e1)])
                try:
                     ov = round(netanalysis.overlap(net, e0, e1), 4)
                except:
                     ov = 0.0
                line = [str(e0), str(e1), str(ov)]
                f.write(' '.join(line) + "\n")
                row = r.readline()
            r.close()
        else:
            for edge in net.edges:
                e0, e1, _ = edge
                e0, e1 = np.sort([e0, e1])
                ov = round(netanalysis.overlap(net, e0, e1), 4)
                line = [str(e0), str(e1), str(ov)]
                f.write(' '.join(line) + "\n")
        f.close()
    else:
        ov_net = pynet.SymmNet()
        for edge in net.edges:
            e_0 = edge[0]
            e_1 = edge[1]
            o = netanalysis.overlap(net, e_0, e_1)
            if o > 0:
                ov_net[e_0, e_1] = o
        return ov_net
    return True


def dics_to_lists(dic_list):
    a = set(key for dic in dic_list for key in dic)
    lists = [[] for _ in dic_list]
    n = len(dic_list)
    for key in a:
        if all(key in dic for dic in dic_list):
            [lists[i].append(dic_list[i][key]) for i in range(n)]
    return lists


def betweeness_centrality(net, edgeBC=True):
    return netext.getBetweennessCentrality(net, edgeBC=edgeBC)

def file_perc(file_path, size=1000):
    giants = []
    sus = []
    nedges = int(list(os.popen("wc -l " + file_path))[0].split()[0])
    e = imap(lineToEdge, open(file_path, 'r'))
    ee = percolator.EvaluationList(e, 'fraclinear', [0.0, 1.0, size], listlen = nedges)
    p = percolator.Percolator(ee, buildNet = False, returnKtree = True)
    for cs in p:
        giants.append(cs.ktree.giantSize)
        sus.append(cs.ktree.getSusceptibility())
    return giants, sus

def plot_perc(giants, label='', line='-', fig=None, ax=None, xlabel=r'$f$', ylabel=r'$P_{GC}$', title='Percolation'):
    size = len(giants)
    nnodes = max(giants)
    f = map(lambda x: x/float(size), range(size))
    n = lambda l: map(lambda x: x/float(nnodes), l)
    if not fig:
        fig, ax = plt.subplots(1)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_title(title)

    ax.plot(f, n(giants), line, label=label)
    ax.legend(loc=0)
    return fig, ax

def weights_to_dic(net):
    dic = {}
    for edge in net.edges:
        e0, e1, w = edge
        e0, e1 = np.sort([e0, e1])
        dic[(e0, e1)] = w
    return dic

# CALL
# giants, sus = applySave(file_perc, ["Edges name"], "output_name.txt")
# f = map(lambda x: x/1000., range(1000))
# n = lambda l : map(lambda x: x/float(nnodes), l)
# plt.plot(f, n(giants), 'r-', f, n(giants_in), 'r--' linewidth=2.0)

############################
# Opverlap
###########################

def get_weight_overlap(net):
    overlap = {}
    for edge in net.edges:
        n1, n2, _ = edge
        n1, n2 = np.sort([n1, n2])
        overlap[n1, n2] = netanalysis.overlap(net, n1, n2)
    return overlap

###3 for printing: obtain overlap according to a net
# Use this with net of residual times, and overlap a dic of
# overlap and total number of calls
def overlap_list(overlap, net, ex=False, all_net=False):
    # extra: return a list of edges as well
    overlaps = []
    weights = []
    extra = []
    for edge in net.edges:
        n1, n2, w = edge
        n1, n2, = np.sort([n1, n2])
        try:
            overlaps.append(overlap[n1, n2])
        except:
            if all_net:
                overlaps.append(np.nan)
        weights.append(w)
        extra.append((n1, n2))
    if ex:
        return overlaps, weights, extra
    else:
        return overlaps, weights

def get_metadata(metadata_path, dic_path, column='age'):
    header = 'uid;age;gender;zipo;lato;longo;subscr;cid;activ_date;disc_date;idc;zipc;latc;longc;adminc;subadminc;localc;distc'
    data = pd.read_csv(metadata_path, header=None, sep=";")
    data.columns = header.split(';')

    dic = pd.read_pickle(dic_path)
    #dic = dict([[v, k] for k, v in dic.items()])

    data['nid'] = data['uid'].apply(lambda x: dic.get(x, -1))
    data = data[data['nid'] != -1]

    return dict(zip(data['nid'], data[column]))

def get_age_difference(edges, dic_path='../data/mobile_network/canarias/canarias_reorder.p', metadata_path='../data/mobile_network/set2/actives_cityloc_200701_200704.dat'):
    dic = get_metadata(metadata_path, dic_path, column='age')
    age_diffs = {}
    for edge in edges:
        e0, e1, w = edge
        try:
            a0, a1 = dic[e0], dic[e1]
            if (a0 > 0) and (a1 > 0):
                age_diffs[e0, e1] = np.abs(a0 - a1)/float(a0 + a1)
        except:
            pass

    return age_diffs


