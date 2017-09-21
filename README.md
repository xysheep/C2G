# Cluster-to-Gate(C2G)
Here is an example of how to run C2G: https://xysheep.github.io/C2G/

Download newest version of C2G: https://github.com/xysheep/C2G/releases

### Quick Start
Given you have one flow or mass cytometry data in fcs format and you want to visualize a kmeans defined cluster, you can do it by following steps.
```MATLAB
[~,markers,~,~,compensated_data] = readfcs_v2(filename);
% Transformation of the data. Cofactor should be 5 for mass cytometry data and 100 for flow cytometry.
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

### Installation


### Hyperparameters

