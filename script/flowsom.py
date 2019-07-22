import FlowCytometryTools
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
import random

from cluster import *
from collections import Counter
from FlowCytometryTools import test_data_dir, test_data_file
from FlowCytometryTools import FCMeasurement
from matplotlib.gridspec import GridSpec
from minisom import MiniSom
from sklearn.cluster import AgglomerativeClustering


class flowsom():
    """"
    Class of flowSOM clustering
    """
    import pandas as pd
    import numpy as np
    import FlowCytometryTools
    from FlowCytometryTools import test_data_dir, test_data_file
    from FlowCytometryTools import FCMeasurement
    from sklearn.cluster import AgglomerativeClustering
    
    def __init__(self, file_address,if_fcs=True, if_drop=True, drop_col=['Time']):
        """
        Read the fcs file as pd.Dataframe

        Parameters
        ----------
        file_address : string 
                       e.g. r'#40 Ab.fcs' or 'flowmetry.csv'
        if_fcs : bool
                 whethe the imput file is fcs file. If not, it should be a csv file
        if_drop : bool
                  define whether some columns should be ignored
        drop_col : list of strings
                   list of column names to be dropped
        """
        if if_fcs:
            self.info = FCMeasurement(ID='Train', datafile=file_address)
            df = self.info.data
        else:
            df = pd.read_csv(file_address)

        self.df = df
        if if_drop:
            self.df = df.drop(drop_col, axis =1)

    def tf(self, tf_str=None, if_fcs=True):
        """
        transform the data, available transform methods include: 'hlog', 'hlog_inv', 'glog', 'tlog'...
        for details, check: https://github.com/eyurtsev/FlowCytometryTools/blob/master/FlowCytometryTools/core/transforms.py#L242

        Parameters
        ----------
        tf_str : string
                 e.g. 'hlog', the transform algorithm
        if_fcs : bool
                 whethe the imput file is fcs file. If not, it should be a csv file
                 only the fcs file could be transformed
                 if it is a csv file, you have to make your own transform function                 
        """
        log_data = pd.DataFrame()
        if if_fcs and tf_str:
            for col in self.df.columns:
                hsample = self.info.transform( tf_str, channels=[col]) # perform transforming for each column
                h_data = hsample.data[col] # get the transforming data
                log_data[col] = h_data.data # store the data into new df
            self.tf_df = log_data
            self.tf_matrix = log_data.values
        else:
            #################################
            ## to-do: transform for ndarray##
            #################################

            self.tf_df = self.df
            self.tf_matrix = self.df.values

    def som_mapping(self, x_n, y_n, d, sigma, lr,batch_size, 
                    neighborhood='gaussian', 
                    tf_str=None, if_fcs=True,
                    seed=10):
        """
        Perform SOM on transform data

        Parameters
        ----------
        x_n : int
              the dimension of expected map
        y_n : int
              the dimension of expected map
        d : int
            vector length of input df
        sigma : float
               the standard deviation of initialized weights
        lr : float 
            learning rate
        batch_size : int
                     iteration times
        neighborhood : string
                       e.g. 'gaussian', the initialized weights' distribution
        tf_str : string
                 tranform parameters, go check self.tf()
                 e.g. None, 'hlog' - the transform algorithm
        if_fcs : bool
                 tranform parameters, go check self.tf()
                 whethe the imput file is fcs file. If not, it should be a csv file
                 only the fcs file could be transformed
                 if it is a csv file, you have to make your own transform function   
        seed : int
               for reproducing
        """
        from minisom import MiniSom

        self.tf(tf_str, if_fcs)
        som = MiniSom(x_n, y_n, d, sigma, lr, neighborhood_function=neighborhood, random_seed=seed) # initialize the map
        som.pca_weights_init(self.tf_matrix) # initialize the weights
        print("Training...")
        som.train_batch(self.tf_matrix, batch_size, verbose=True)  # random training
        print("\n...ready!")
        self.x_n = x_n
        self.y_n = y_n
        self.map_som = som
        self.weights = som.get_weights()
        self.flatten_weights = self.weights.reshape(x_n*y_n, d)
        
    def meta_clustering(self, cluster_class, min_n, max_n, iter_n, 
                        resample_proportion=0.7,verbose=False):
        """
        Perform meta clustering on SOM

        Parameters
        ----------
        cluster_class : class
                        e.g. KMeans, a cluster class, like "from sklearn.cluster import KMeans"
        min_n : int
                the min proposed number of cluster
        max_n : int
                the max proposed number of cluster
        iter_n : int
                 the iteration times for each number of clusters
        resample_proportion : float
                              within (0, 1), the proportion of re-sampling when computing clustering
        verbose : bool
                  whether print out the clustering process
        """
        
        # initialize cluster
        cluster_ = ConsensusCluster(cluster_class, 
                                    min_n, max_n, iter_n, 
                                    resample_proportion=resample_proportion)
        cluster_.fit(self.flatten_weights, verbose) # fitting SOM weights into clustering algorithm

        self.cluster_map = cluster_
        self.bestk = cluster_.bestK # the best number of clusters in range(min_n, max_n)

        # get the prediction of each weight vector on meta clusters (on bestK)
        self.flatten_class = cluster_.predict_data(self.flatten_weights)
        self.map_class = self.flatten_class.reshape(self.x_n, self.y_n) 

    def vis(self, t, with_labels, node_size, edge_color):
        """
        Visualize the meta cluster result with minimal spanning tree

        Parameters
        ----------
        t : int
            total number of nodes, n = t * bestK
        with_labels : bool
                      whether the node will be visualized with its cluster label
        node_size : int
        edge_color : string
                     e.g 'b', the color of edges
        """
        from matplotlib.gridspec import GridSpec
        import networkx as nx
        import numpy as np
        from collections import Counter
        import matplotlib.pyplot as plt
        import matplotlib.cm as cm
        
        # generate n clusters (n = bestK * t)
        self.cluster_map.bestK = self.bestk * t
        self.over_class = self.cluster_map.predict_data(self.flatten_weights)
        centroid_list = []
        
        # Compute the centroid for each clusters
        for i in np.unique(self.over_class):
            centroid = np.mean( self.flatten_weights[self.over_class == i], axis=0 )
            centroid_list.append(centroid)
        self.centroid_list = centroid_list

        # Generate a fully connected graph of n cluster centroids for future computation 
        # on minimal spanning tree
        # (node: centroid of cluster, weight of edge: distance between two nodes)
        G = nx.Graph()
        
        for i in range(len(centroid_list)):
            for j in range(i+1, len(centroid_list)):
                # compute the distance between two nodes
                w = np.sqrt(np.dot(centroid_list[i], centroid_list[i]) - 
                            2 * np.dot(centroid_list[i], centroid_list[j]) + 
                            np.dot(centroid_list[j], centroid_list[j]))
                w /= 1
                G.add_edge(i, j, weight=w)
        self.graph = G
        mst = nx.minimum_spanning_tree(G) # compute the minimal spanning tree graph
        self.mst = mst
        
        # generate the plot
        edges,weights = zip(*nx.get_edge_attributes(mst,'weight').items())
        print (self.bestk)
        color_list = cm.rainbow(np.linspace(0, 1, self.bestk))
        color_map = []
        for node in mst:
            class_id, _ = Counter(self.flatten_class[self.over_class == node]).most_common()[0]
            try:
                color_map.append(color_list[class_id])
            except:
                print ('something wrong with plotting cluster %d!' % class_id)
        nx.draw(mst, 
                with_labels=with_labels, 
                node_size = node_size,
                node_color=color_map, 
                edgelist=edges, 
                edge_color=edge_color, 
                width=weights*100, 
                edge_cmap=plt.cm.Blues)
        # splt.show()
  
    def labeling(self, verbose=True):
        """
        Make prediction for the whole dataset, add a column of 'category' with prediction. 
        """
        label_list = []
        for i in range(len(self.tf_matrix)):
            # print the milestone
            if verbose:
                if i % 10000 == 0:
                    print ('%d samples done...' % i)

            xx = self.tf_matrix[i, :] # fetch the sample data
            winner = self.map_som.winner(xx) # make prediction, prediction = the closest entry location in the SOM
            c = self.map_class[winner] # from the location info get cluster info
            label_list.append(c)
        self.df['category'] = label_list
        self.tf_df['category'] = label_list