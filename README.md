# FlowSOM
This repository contains a Python implementation of [FlowSOM](http://bioconductor.org/packages/release/bioc/html/FlowSOM.html) algorithm for clustering and visualizing a mass cytometry data set. 

# Installation
Download FlowSOM to a directory of your choice and then run:

    pip install -r requirements.txt
    
    
How to use it
------------------
In order to use FlowSOM you need your data saved as a .csv file or a .fcs file.
```python
file = r'flowmetry.fcs'
```
or 
```python
file = 'flowmetry.csv'
```

Then you can run FlowSOM just as follows:
```python
from flowSOM import *
fsom = flowsom(file) # read the data
fsom.som_mapping(50, 50, 31, sigma=2.5, 
                 learning_rate=0.1, batch_size=100)  # trains SOM with 100 iterations
fsom.meta_clustering(AgglomerativeClustering, min_n=40, 
                     max_n=45, 
                     iter_n=3) # train the meta clustering for cluster in range(40,45)       
```

### Use the trained output

After the training, you will be able to:

* Get the weights of SOM with method `fsom.map_som`
* Get the best number of clustering with method `fsom.bestk`
* Get the prediction dataframe with method `fsom.df` and `fsom.tf_df`
* Visualize the final clustering outcome with method`fsom.vis`

Examples
-------------------------
The demo code could be found [here](https://github.com/Hatchin/FlowSOM/blob/master/demo.ipynb).

The distance map of SOM trained from a sample flow cytometry [data](https://github.com/Hatchin/FlowSOM/blob/master/flowmetry_transformed.csv):

<img src="https://github.com/Hatchin/FlowSOM/blob/master/som.png" alt="Flow example">

The visualization example after meta-clustering using Minimal Spanning Tree (MST):
<img src="https://github.com/Hatchin/FlowSOM/blob/master/mst.png" alt="MST example">

FlowSOM Algorithm
--------------------------

FlowSOM analyzes flow or mass cytometry data using a self-Organizing Map (SOM). Using a two-level clustering and star charts, FlowSOM helps to obtain a clear overview of how all markers are behaving on all cells, and to detect subsets that might be missed otherwise. 

The algorithm consists of four steps: 
    - reading the data;
    - building a Self-Organizing Map;
    - building a minimal spanning tree;
    - computing a meta-clustering. 
    
### Self-Organizing Map
SOM is a type of unsupervised Artificial Neural Network able to convert complex, nonlinear statistical relationships between high-dimensional data items into simple geometric relationships on a low-dimensional display. [Introduction](https://heartbeat.fritz.ai/introduction-to-self-organizing-maps-soms-98e88b568f5d)

### Minimum Spanning Tree
```
A minimum spanning tree (MST) or minimum weight spanning tree is a subset of the edges of a connected, edge-weighted undirected graph that connects all the vertices together, without any cycles and with the minimum possible total edge weight.
```

### Meta-clustering
The meta-clustering technique conducted on the SOM is hierarchical consensus meta-clustering, which clusters the weights of trained SOM into different groups. 

Acknowledge
-----------------
FlowSOM is built based on [FlowCytometryTools](https://github.com/eyurtsev/FlowCytometryTools), [MiniSom](https://github.com/JustGlowing/minisom) and [Consensus Clustering](https://github.com/ZigaSajovic/Consensus_Clustering).
