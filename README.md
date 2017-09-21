# Cluster-to-Gate(C2G)
Here is an example of how to run C2G: https://xysheep.github.io/C2G/

Download newest version of C2G: https://github.com/xysheep/C2G/releases

### Quick Start
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

### Installation
In general, you need to have MATLAB installed and directory of "libs" and "src" in your MATLAB path. You can ensure it by following steps. 
1. Have a MATLAB installed on your computer
2. Download the most recent release from [here](https://github.com/xysheep/C2G/releases)
3. Unzip the downloaded folder to the place you like
4. In MATLAB, go to the unzipped folder
5. Right click the "libs" and "src" folder and choose "add to path"
### Hyperparameters

