import numpy as np
from analytics import *
import matplotlib.pyplot as plt; plt.ion()


class TieAnalytics(object):
    def __init__(self, nodes, times=None, query='', graph):
        """
        times - dict with
        """
        self.nodes = nodes
        if times is None:
            assert query, "times or query must be provided"
            self.times = self._query_to_times(query)
        else:
            self.times = times
        self.stats = {}

    def _get_active_times(self):
        pass


    def _query_to_times(self, query):
        pass


    def get_iet_distribution(self):
        pass


    def get_overlap(self):
        pass


    def get_bursty_cascades(self, dt=365):
        pass


    def get_reciprocity(self):
        pass

