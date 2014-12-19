Author: Zhiqiang Xu, 05/2012
Email : zxu1@e.ntu.edu.sg

This package contains the MATLAB implementation of the Bayesian Attributed Graph Clustering (BAGC) algorithm, proposed in the paper "A Model-based Approach to Attributed Graph Clustering" (Xu et al., 2012).

Note that this code is only applicable to unweighted and undirected graphs with vertices associated with categorical attributes. The data format for adjacency matrix is coordinate list (COO) and for attributes a table with a row for each vertex and a column for each attribute. Note that Matlab array indices are "1-based". In addition to data, the input also includes the cluster number, the number of optimizing iterations which is set to 5 by default, and the initial clustering (partition, or assignment) of vertices optionally provided by users. The user-provided initial clustering may mitigate the cost of time and memory caused by the default initialization which is the same as in competitor used in paper, Inc-cluster, for fairness. The output includes modularity scalar, entropy vector with an element for each attribute, and the optimizing time.

The following files are included:

1. An example attributed graph with 100 vertices
- adj_coo.txt    :    a 100x100 adjacency matrix stored in COO (coordinate list) format.
- attr_tab.txt   :    a 100x2 attribute table with a row for each vertex and a column for each attribute.

2. MATLAB codes
- example.m      :    an example that illustrates the use of the code on the given attributed graph. The output clustering has a modularity value 0.5513, and entropies 0.2318 for attribute 1 (column 1 of attr_tab.txt) and 0.3849 for attribute 2 (column 2 of attr_tab.txt). The running time is about 0.005s on a machine with 32-bit Windows OS, 2.67GHz 4-core CPU and 4G memory.
- bagc.m         :    the main function for BAGC.
- formatted.m    :    converts raw data to an internal representation.
- init.m         :    initializes the posterior distribution of the cluster labels.
- optimize.m     :    iteratively optimizes the posterior of the cluster labels.
- evaluate.m     :    evaluates the structure and attribute quality of the resulting clustering.
- run.m		 :    provide a single function for the use of the code.

3. Documentation
- readme.txt     :    this file.