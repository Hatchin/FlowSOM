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

