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

Examples
-------------------------
The demo code could be found [here](https://github.com/Hatchin/FlowSOM/blob/master/demo.ipynb).

The distance map of SOM trained from a sample flow cytometry [data](https://github.com/Hatchin/FlowSOM/blob/master/flowmetry_transformed.csv):

<img src="https://github.com/Hatchin/FlowSOM/blob/master/som.png" alt="Flow example">

The visualization after meta-clustering using Minimal Spanning Tree (MST):
<img src="https://github.com/Hatchin/FlowSOM/blob/master/mst.png" alt="MST example">


