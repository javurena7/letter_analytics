# letter_analytics
**Network analysis tools for communication logs**

- *analytics.py*: basic analytics functions used by one (or more) views. 
- *plots.py*: basic plotting functions for other cases [in development]

Main perspectives on how to analyze networks.
 - *ego.py*: main classes for anlysis of ego networks. Including basic ego statistics and social signatures. 
 - *network_stats.py*: main classes for network analysis [in development]
 - *place.py*: main classes for analyzing geographical data [in development]
 - *tie.py*: main classes for analyzing tie-level data [in development]

## Example
Basic use of *ego.py*, where `query` is a json query result for an ego network including incoming and outgoing letters to/from `node`, including timestamps.

The two main temporal analytics of the data are the social signature (which tells how people divide their communication patterns) and network attractiveness, which includes basic measures of how much a node interacts with the network vs how much people in the network interact with the node. 

### Social Signature
```python
ea = ego.EgoAnalytics(node, query=query)
ea.social_signature(bin_type='year', bin_n=5, max_rank=20) #obtain signatures for 5-year bins, where signatures have at most 20 people
ea.ss.t_signatures #dict containing the signatures and the nodes which appear in each
ea.ss.average_signature #list of average signature
ea.ss.plot_all_signatures()

ea.ss.get_year_edges() #Obtain the year edges used in binning

ea.ss.update(bin_type='linear', bin_n=10) #Get new signatures where the timespan of interactions is linearly split into 10 bins
```
### Network Attractiveness
```python
ea.dynamic_attractiveness(bin_type='year', bin_n=5)
ea.da.t_stats #Obtain the dict of stats per bin
ea.da.update(bin_type='linear', bin_n=10)
```
