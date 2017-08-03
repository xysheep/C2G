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
m.view_gates(data,markernames,'n_lines',1);
outtable = m.show_f_score(label); 

%% Load CYTOF datasets
% Read multiple fcs files. Each fcs file correspond to one cell population
% Automatically generate cell labels based on fsc files. FCS files are
% store in the "test" folder. "CD4_Effmen.fcs", "CD4_naive.fcs", and
% "CD8_naive.fcs" are cells of target populations and "ctr.fcs" contain all
% cells. Only work with surface protein markers. 
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
% Precluster the ungated cells

rng(9464)
%label = cluster_ungated_gmm(data,ori_l,100,60);
label = cluster_ungated(data,ori_l);
% Perform the anlysis
m_ori = C2G(data,label,ori_l,'markernames',markers);
%m = C2G2(data,label,ori_l,means,covs,overlap2d);
% Record the fareast cell from the normal distribution. 
%m.draw_gates(data,label, density);
m_ori.view_gates(data,markers,'n_lines',3,'ignore_small',0);
m_ori.show_f_score(ori_l);

%% Consistentcy to subsampling
n_sub = 10;
sub_portion = 0.9;
n_event = size(data,1);
m_sub = cell(n_sub,1);
nmi = cell(n_sub,1);
rng(9464)
for i = 1:n_sub
    sub_idx = randsample(n_event,round(n_event*0.9));
    tmp_d = data(sub_idx,:);
    tmp_label = ori_l(sub_idx);
    label = cluster_ungated(tmp_d,tmp_label);
    m_sub{i} = C2G(tmp_d, label, tmp_label,'markernames',markers);
    
    [~,nmi{i}] = m_sub{i}.show_f_score(tmp_label);
    drawnow;
end
for i = 1:n_sub
    m_sub{i}.view_gates(tmp_d,markers,'n_lines',1,'onepanel',true);
end
%% Generate gating hierarchy for K-means defined populations (K=10)
rng(9464)
km_l = kmeans(data,10);
% Precluster the ungated cells
label = cluster_ungated(data,km_l);
new_km_l = km_l;
% Perform the anlysis
m_km = C2G(data,label,new_km_l,'markernames',markers);
m_km.view_gates(data,markers,'n_lines',5,'ignore_small',0);
m_km.show_f_score(new_km_l);

%% Visualization of SPADE overclustered results
[clustered_data,marker]  = readfcs_v2('testdata/spade/exported_results/clustering_result_added_in_fcs/ctr.fcs');
surface_idx = [3 4 6 8 9 11 12 13 22 24 25 27];
ori_data = clustered_data(surface_idx,:)';
data =  flow_arcsinh(ori_data,5);
label = clustered_data(end,:)';
tic;m_spade = C2G(data, label, label,'markernames',marker(surface_idx));toc;
%m.view_gates(d,marker(surface_idx),'ignore_small',10);
m_spade.show_f_score(label);



%% Scaffold CYTOF data
% 21 protein markers
[ori_data, marker] = readfcs_v2('testdata/bigdata/TIN_BLD1_Untreated_Day3.fcs');
data = flow_arcsinh(ori_data,5);
marker_idx = [10,16,17,18,19,21,22,24,29,31,32,33,39,40,41,46,47,49,50,51,52];
d = data(marker_idx,:)';
rng(9464);
label = kmeans(d,10);
tic;m = C2G(d, label, label,'markernames',marker(marker_idx));toc;
m.view_gates(d,marker(marker_idx),'ignore_small',3000,'n_lines',4);
m.show_f_score(label);