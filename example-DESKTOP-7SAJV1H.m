%% Generating 2D Gating Hierarhcy from Clustered Cytometry Data
% Our algorithm cluster-to-gates(C2G) can visualize a clustered cytometry
% data in 2D gating hierarchy. The overall input to C2G are two variables:
% "data" and "label". "data" is a N-by-M matrix where N is the number of
% cells and M is the number of markers. "label" is a N-by-1 matrix and this
% matrix represent which cluster each cell belongs to (0 for ungated
% cells).  In this page, we will demonstrate how to use C2G to generate
% gating hierarchy on a simulated data and on a CYTOF dataset that are
% manually gated.
% The software and testing data can be downloaded *here* . 

%% Initialization and load the simulated data
% This simulated data set contain 20000 cells, which consist of 5
% populations.Two of the populations are addressed as known and labeled as
% 1 and 2. The cells in other populations are labeled 0 which means
% "ungated". A visualization of the original data is also shown in the
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
% have no "ungated" cells, preclustered step can be skipped. If you skip
% precluster, the second and third parameter in function "C2G" should be
% the same. Compute local density step is also optional, it can decrease
% time cost in computing density when you want to try different parameter
% setting in "C2G". If you skipped it, no need to input the fourth
% parameter in "C2G".

% Precluster the simulated data
rng(9464);
preclustered_label = cluster_ungated(data,label);
figure('Position',[680 478 560 420]);
col = [0.8 0.8 0.8;hsv(length(unique(preclustered_label)))];
scatter3(data(:,1),data(:,2),data(:,3),1,col(preclustered_label+1,:));
xlabel('Marker 1')
ylabel('Marker 2')
zlabel('Marker 3')
axis([-5 15 -5 15 -5 15]);

% Call main part of the program and return a object m that store the
% results. The option "ignore_ratio" mean percentage of low density cells
% ignored when compute overlap between different populations

m = C2G(data,preclustered_label,label,'markernames',markernames,'color',col); 
%m = C2G2(data,preclustered_label,label,means,covs,overlap2d,'markernames',markernames,'color',col); 
% Draw the obtained gating hierarchy
% Show statistics
%m.draw_gates(data,preclustered_label, density);
m.view_gates(data,markernames,20);
outtable = m.show_f_score(label); 

%% Load CYTOF datasets
% Read multiple fcs files. Each fcs file correspond to one cell population
% Automatically generate cell labels based on fsc files. FCS files are
% store in the "test" folder. "CD4_Effmen.fcs", "CD4_naive.fcs", and
% "CD8_naive.fcs" are cells of target populations and "ctr.fcs" contain all
% cells. Only work with surface protein markers. 
fdname = 'testdata';
[ori_data,ori_l,ori_markers]=load_mul_fcs(fdname,'ctr.fcs');

surface_idx = [3 4 6 8 9 11 12 13 22 24 25 27];
data = ori_data(:,surface_idx);
markers = ori_markers(surface_idx);
n_markers = length(markers);

%% Generate gating hierarchy for manually gated populations
% Precluster the ungated cells
rng(9464)
%label = cluster_ungated_gmm(data,ori_l,100,60);
label = cluster_ungated(data,ori_l);
% Perform the anlysis
m = C2G(data,label,ori_l,'markernames',markers);
%m = C2G2(data,label,ori_l,means,covs,overlap2d);
% Record the fareast cell from the normal distribution. 
%m.draw_gates(data,label, density);
m.view_gates(data,markers,20);
m.show_f_score(ori_l);

%% Generate gating hierarchy for K-means defined populations (K=10)
rng(9464)
km_l = kmeans(data,10);
% Precluster the ungated cells
label = cluster_ungated(data,km_l);
new_km_l = km_l;
% Perform the anlysis
m = C2G(data,label,new_km_l,'markernames',markers);
m.view_gates(data,markers);
m.show_f_score(new_km_l);



%% Perform on large CYTOF data
% This dataset has 200k cells and 10 protein markers
load testdata\Screen_PACA-1-3_008.mat
ori_data = data;
data = ori_data(7:16,ori_data(21,:)>0.1)';
markers = marker_names(7:16);

rng(9464)
km_l = kmeans(data,10);
% Precluster the ungated cells
label = cluster_ungated(data,km_l);
new_km_l = km_l;
% Perform the anlysis
m = C2G(data,label,new_km_l,'ignore_ratio',0.2);
m.view_gates(data,markers);
m.show_f_score(new_km_l);

% Pretent not know some clusters
km_l(km_l > 6) = 0;
label = cluster_ungated(data,km_l);
new_km_l = km_l;
% Perform the anlysis
m = C2G(data,label,new_km_l,'ignore_ratio',0.2);
m.view_gates(data,markers);
m.show_f_score(new_km_l);
% Around 10 minutes, f-score 0.86.
% bottlenect is convert cell into grid.

%% Scaffold CYTOF data
% 21 protein markers
[ori_data, marker] = readfcs_v2('testdata/TIN_BLD1_Untreated_Day3.fcs');
data = flow_arcsinh(ori_data,5);
marker_idx = [10,16,17,18,19,21,22,24,29,31,32,33,39,40,41,46,47,49,50,51,52];
d = data(marker_idx,:)';
rng(9464);
label = kmeans(d,10);
tic;m = C2G(d, label, label,'markernames',marker(marker_idx));toc;
m.view_gates(d,marker(marker_idx),'ignore_small',1000);
m.show_f_score(label);