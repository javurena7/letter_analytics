import numpy as np
import pandas as pd
import networkx as nx
from analytics import *
import matplotlib.pyplot as plt; plt.ion()

logs = '1500_1990/30/times_dic.txt'

class NetAnalytics(object):
    def __init__(self, net):
        self.net = net

    def net_attractiveness(self):
        """
        Obtain measures of network attractiveness from Wang 2016. 
        
        Returns:
        stats: dictionary where keys are node ids, values are dicts
        """
        stats = {}
        for node in self.net:

            in_set = {k: self.net.get_edge_data(k, node)['weight'] for k in net.predecessors(node)}
            out_set = {k: self.net.get_edge_data(node, k)['weight'] for k in net.successors(node)}
            stats[node] = ego_reciprocity(in_set, out_set)

        return stats


    def percolation(self, p_list):
        pass

