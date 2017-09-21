# Cluster-to-Gate(C2G)
Cluster-to-Gate is a tool developed to visualize clustering results of multi-dimensional flow/mass cytometry data. 

Here is an example of how to run C2G: https://xysheep.github.io/C2G/. 

Download newest version of C2G: https://github.com/xysheep/C2G/releases

## Quick Start
If you already have MATLAB installed on your computer you can try the quick start. Given you have one flow or mass cytometry data in fcs format and you want to visualize cell populations defined by kmeans, you can do it using the following code. In the transformation step, cofactor should be 5 for mass cytometry data and 100 for flow cytometry data.
```MATLAB
[~,markers,~,~,compensated_data] = readfcs_v2(filename);
% Transformation of the data. 
transformed_data = flow_arcsinh(compensated_data, cofactor)';
% Run k-means where k equal to 10
label = kmeans(transformed_data, 10);
% The main function of C2G
m = C2G(transformed_data, label, label);
% Visualize the gating hierarchy
m.view_gates(transformed_data, markers);
% See statistics of the results
m.show_f_score(label);
```
Alternatively, if you start from a data matrix already loaded into MATLAB and have a clustering results defined by other programs, you can run C2G in the following way. Here, you data should be M-by-N matrix where M is number cells while N is number markers. The clustering label is stored in a N-by-1 vector. Clustering label can be any non-negative integers and 0 means unlabeled.  
```MATLAB
transformed_data = flow_arcsinh(data, cofactor);
% Pre-cluster the unlabeled cells
newlabel = cluster_ungated(data,label);
% The main function of C2G
m = C2G(transformed_data, newlabel, label);
% Visualize the gating hierarchy
m.view_gates(transformed_data);
% See statistics of the results
m.show_f_score(label);
```
## Installation
In general, you need to have MATLAB installed and directory of "libs" and "src" in your MATLAB path. You can ensure it by following steps. 
1. Have a MATLAB installed on your computer
2. Download the most recent release of C2G from [here](https://github.com/xysheep/C2G/releases)
3. Unzip the downloaded folder to the place you like
4. In MATLAB, go to the unzipped folder
5. Right click the "libs" and "src" folder and choose "add to path"
## Parameters
C2G has the following hyperparameter that you can tune according to the shape of your data. 
### function cluster_ungated
**'pre_cluster_perc'** Minimum percentage of unlabeled cells that should be pre-clustered. Default is 0.95. 

**'outlier_level'** Percentage of low density cells of target populations ignored when decide whether a unlabeled cells overlaped with target populations. Default is 0.05. 
### function C2G
**'showdetail'** Whether to show the f-score of each cell populations after each iteration. Default is true. 
**'trivial_gate'** Smallest number of cells a single gate must exclude when it's the only gate in that marker pair. Default is 50.
### function gatingTree.view_gates
**'fontsize'** Default is 20. 

**'n_lines'** Number of rows in output figures that have multiple panels. Default is 3. 

**'ignore_small'** This is a threshold. When a node has only one gate and this gate fail to separate more than 'ignore_small' number of cells, this gate will not be drawn.  This function is useful to filter out some unnecessary gates and only leave the important gates that show the hierarchy.  Compare to 'trivial_gate' parameter in C2G section, this parameter only affect the visualization. Default is 100. 

**'onepanel'** Whether the gating tree is present in the same panel of the gating sequences. Default is false. 
