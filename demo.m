%% Generating 2D Gating Hierarhcy from Clustered Cytometry Data
% Our algorithm cluster-to-gates(C2G) can visualize a clustered cytometry
% data in 2D gating hierarchy. The overall input to C2G are two variables:
% "data" and "label". "data" is a N-by-M matrix where N is the number of
% cells and M is the number of markers. "label" is a N-by-1 matrix and this
% matrix represent which cluster each cell belongs to (0 for ungated
% cells).  In this page, we will demonstrate how to use C2G to generate
% gating hierarchy on a simulated data and on a CyTOF dataset that are
% manually gated. The software and testing data can be downloaded
% <https://github.com/xysheep/C2G/releases here> .

%% Initialization and load the simulated data
% This simulated data set contains 20000 cells, which consists of 5
% populations.Two of the populations are considered as target and labeled as
% 1 and 2. The cells in other populations are labeled 0 which means
% "unlabeled". A visualization of the original data is also shown in the
% following sections.
addpath('src')
addpath('libs')
load('testdata/simulated.mat','data','label');
markernames ={'Marker 1','Marker 2','Marker 3'};
fprintf('Size of "data" is %d-by-%d\n',size(data));
fprintf('Size of "label" is %d-by-%d\n',size(label));

%% Visualize the simulated data and pre-clustered simulated data in 3D
% In real application, this section is not necessay.
col = [0.8 0.8 0.8;hsv(length(unique(label))-1)];
figure('Position',[680 478 560 420]);
scatter3(data(:,1),data(:,2),data(:,3),1,col(label+1,:));
xlabel('Marker 1')
ylabel('Marker 2')
zlabel('Marker 3')
axis([-5 15 -5 15 -5 15]);

%% Run the analysis on the simulated data
% This section perform the analysis on the simulated data. If your data
% have no "unlabeled" cells, precluster step can be skipped. If you skip
% precluster, the second and third parameter in function "C2G" should be
% the same. 

% Precluster the simulated data
rng(9464);
preclustered_label = cluster_ungated(data,label);

% Call main part of the program and return a object m that store the
% results. 
m = C2G(data,preclustered_label,label,'markernames',markernames); 
% Draw the obtained gating hierarchy
m.view_gates(data,markernames,'n_lines',1,'onepanel',true);
% Show statistics
outtable = m.show_f_score(label); 

%% Load a moderate-sized CyTOF dataset with ~10,000 cells
% This is a dataset about T cell signaling (Krishnaswamy et al, Science,
% 2014). It can be downloaded from
% <https://www.c2b2.columbia.edu/danapeerlab/html/dremi-data.html here>.
% This section will read multiple fcs files. Each fcs file correspond to
% one cell population manually gated by the authors. C2G will automatically
% generate cell labels based on fsc files. FCS files are store in the
% "testdata" folder. "CD4_Effmen.fcs", "CD4_naive.fcs", and "CD8_naive.fcs"
% are cells of target populations (manually defined) and "ctr.fcs" contain
% all cells. In this example, we only use surface protein markers.
clear
close all
addpath('src')
addpath('libs')
fdname = 'testdata';
[ori_data,ori_l,ori_markers]=load_mul_fcs(fdname,'ctr.fcs');

surface_idx = [3 4 6 8 9 11 12 13 22 24 25 27];
data = ori_data(:,surface_idx);
markers = ori_markers(surface_idx);
n_markers = length(markers);

%% Generate gating hierarchy for manually gated populations
% In this section, we'll use the data loaded from previous section and
% treat the three manually defined cell populations as target populations.
% Here, number of "unlabeled" cells is around two-fold of cells in target
% populations.

% Precluster the ungated cells
rng(9464)
label = cluster_ungated(data,ori_l);
% Perform the anlysis
m_ori = C2G(data,label,ori_l,'markernames',markers);
% Visualize the results
m_ori.view_gates(data,markers,'n_lines',3,'ignore_small',0,'onepanel',true);
m_ori.show_f_score(ori_l);

%% Generate gating hierarchy for K-means defined populations (K=10)
% In this section, we will use the same CyTOF data used in previous
% section. The cell labels are defined by k-means algorithm where k equal
% to 10. All of the 10 clusters will be considered as target populations,
% in other words,there're no "unlabeled" cells.

rng(9464)
km_l = kmeans(data,10);
new_km_l = km_l; % Since all populiations is known, no need to pre-cluster
% Perform the anlysis
m_km = C2G(data,km_l,km_l,'markernames',markers);
% Visualize the results
m_km.view_gates(data,markers,'n_lines',4,'ignore_small',300);
m_km.show_f_score(km_l);

%% Apply C2G to large CyTOF dataset with ~140,000 cells
% This is a larger CyTOF dataset with 140k cells and 21 protein markers
% (Spitzer, Matthew H., et al, Science, 2015 ). This data can be download
% free from <https://community.cytobank.org/cytobank/experiments/60068
% here>. (The dataset is hosted on Cytobank and you need to register to
% download it) This dataset is clustered by k-means where k equal to 10.
% All of the 10 k-means defined clusters are treated as target populations.
% This part will take around 10 minutes on a desktop.
[ori_data, marker] = readfcs_v2('testdata/bigdata/TIN_BLD1_Untreated_Day3.fcs');
% Transform the CyTOF
data = flow_arcsinh(ori_data,5);
% Select protein markers
marker_idx = [10,16,17,18,19,21,22,24,29,31,32,33,39,40,41,46,47,49,50,51,52];
d = data(marker_idx,:)';
rng(9464);
label = kmeans(d,10);
% Perform the anlysis
tic;m = C2G(d, label, label,'markernames',marker(marker_idx));toc;
% Visualize the results
w = warning ('off','all');
m.view_gates(d,marker(marker_idx),'ignore_small',3000,'n_lines',4);
warning(w);
m.show_f_score(label);

%% Reference
% Krishnaswamy, Smita, Matthew H. Spitzer, Michael Mingueneau, Sean C. Bendall, Oren Litvin, Erica Stone, Dana Pe�er, and Garry P. Nolan. "Conditional density-based analysis of T cell signaling in single-cell data." _Science_ 346, no. 6213 (2014): 1250689.
%
% Spitzer, Matthew H., Pier Federico Gherardini, Gabriela K. Fragiadakis, Nupur Bhattacharya, Robert T. Yuan, Andrew N. Hotson, Rachel Finck et al. "An interactive reference framework for modeling a dynamic immune system." _Science_ 349, no. 6244 (2015): 1259425.