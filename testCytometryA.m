%% This test is used to generate all figures in manuscript
addpath('src')
addpath('libs')

%% Figure 2 work flow of method
load('testdata/simulated.mat','data','label');
markernames ={'Marker 1','Marker 2','Marker 3'};
fprintf('Size of "data" is %d-by-%d\n',size(data));
fprintf('Size of "label" is %d-by-%d\n',size(label));

preclustered_label = cluster_ungated(data,label);
col = [0.8 0.8 0.8;hsv(length(unique(preclustered_label))-1)];
% print('manuscript/Figures/figure 2b.pdf','-dpdf');

downidx = downsample(1:size(data,1),20);
fa=figure('Position',[680 678 280 210]);
scatter3(data(downidx,1),data(downidx,2),data(downidx,3),1,col(label(downidx)+1,:));
xlabel('Marker 1')
ylabel('Marker 2')
zlabel('Marker 3')
axis([-5 15 -5 15 -5 15]);
title(sprintf('Original Data \nwith 2 Target Population'))
% print('manuscript/Figures/figure 2a.pdf','-dpdf');

fc=figure('Position',[680 678 280 210]);  % This figure might be drawn in each marker separately
scatter3(data(downidx,1),data(downidx,2),data(downidx,3),1,col(preclustered_label(downidx)+1,:));
xlabel('Marker 1')
ylabel('Marker 2')
zlabel('Marker 3')
axis([-5 15 -5 15 -5 15]);
title('Preclustered Data')
% print('manuscript/Figures/figure 2c.pdf','-dpdf');
% Subpopulation in first gate
downidx12 = downidx(ismember(label(downidx),[1 2]));
fh=figure('Position',[680 678 280 210]);  % This figure might be drawn in each marker separately
scatter3(data(downidx12,1),data(downidx12,2),...
    data(downidx12,3),1,col(preclustered_label(downidx12)+1,:));
xlabel('Marker 1')
ylabel('Marker 2')
zlabel('Marker 3')
axis([-5 15 -5 15 -5 15]);
title('Preclustered Data')
% print('manuscript/Figures/figure 2h.pdf','-dpdf');
% Main C2G
[m,fdemo] = C2G(data,preclustered_label,label,'color',col,'markernames',markernames);
% print(fdemo{1},'manuscript/Figures/figure 2d.pdf','-dpdf')
% print(fdemo{2},'manuscript/Figures/figure 2i.pdf','-dpdf')
m.view_gates(data,markernames,'n_lines',1,'onepanel',true);
% Show statistics
outtable = m.show_f_score(label); 

%% Load Smaller CyTOF data
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

%% Figure 3, Table 1, Manual Gating
rng(9464)
label = cluster_ungated(data,ori_l);
% Perform the anlysis
m_ori = C2G(data,label,ori_l,'markernames',markers);
% Visualize the results
m_ori.view_gates(data,markers,'n_lines',3,'ignore_small',0,'onepanel',true);
m_ori.view_gates_contour(data,markers,'n_lines',2,'ignore_small',0,'onepanel',true);
m_ori.show_f_score(ori_l);
% manu_flowtype = compareflowtype(data, label);

%% Figure 4,  Table S1, kmeans
rng(9464)
km_l = kmeans(data,10);
new_km_l = km_l; % Since all populiations is known, no need to pre-cluster
% Perform the anlysis
m_km = C2G(data,km_l,km_l,'markernames',markers,'trivial_gate',200,'ratio_trivial_gate',0.1);
% Visualize the results
m_km.view_gates_contour(data,markers,'n_lines',4,'ignore_small',0);
m_km.show_f_score(km_l);
%kmeans_flowtype = compareflowtype(data, new_km_l);
%% Figure 5, Table 2, SPADE
d = readfcs_v2('testdata\SPADE\exported_results\clustering_result_added_in_fcs\spade1.fcs');
SPADE_l = d(end,:)';

m_spd = C2G(data, SPADE_l, SPADE_l, 'markernames',markers, 'trivial_gate', 100, 'grid_size', 30);
m_spd.view_gates(data,markers,'ignore_small',0,'onepanel',true);
m_spd.show_f_score(SPADE_l);


%% Figure 6, knowone vs. knowall
rng(9464)
km_l = kmeans(data,10);
% Know one case
out_knowone=cell(10,1);
parfor i = 1:10
    new_km_l = km_l;
    new_km_l(~ismember(new_km_l,i)) = 0;
    label = cluster_ungated(data,new_km_l,'cluster_amount',0.9);
    m=C2G(data,label,new_km_l);
    out_knowone{i}=m.show_f_score(new_km_l);
end
% Know all case
rng(9464)
label = cluster_ungated(data,km_l,'cluster_amount',0.9);
new_km_l = km_l;
m=C2G(data,label,new_km_l);
out_knowall = m.show_f_score(new_km_l);
% Know 3 case
out_knowthree = cell(10,1);
rng(42)
parfor i = 1:10
    know_idx = [i-1 mod(i,10) mod(i+1,10)]+1;
    new_km_l = km_l;
    new_km_l(~ismember(new_km_l,know_idx)) = 0;
    label = cluster_ungated(data,new_km_l,'cluster_amount',0.9);
    m = C2G(data,label,new_km_l);
    out_knowthree{i} = m.show_f_score(new_km_l);
end
tmp_know1 = cell2mat(out_knowone);
overall_score = zeros(10,3);
tmp_know3 = cell2mat(out_knowthree);
save know1310.mat
for i = 1:10
    overall_score(i,1) = tmp_know1(i,6);
    overall_score(i,2) = mean(tmp_know3(tmp_know3(:,1)==i,6));
    overall_score(i,3) = out_knowall(i,6);
end
figure
hb = bar(overall_score);
set(hb(1), 'FaceColor',[0.25 0.25 0.25], 'EdgeColor', [0 0 0])
set(hb(2), 'FaceColor',[1 1 1], 'EdgeColor', [0 0 0])
set(hb(3), 'FaceColor',[0.75 0.75 0.75], 'EdgeColor', [0 0 0])
axis([0 11 0 1.15])
legend({'Know one','Know three','Know all'},'Location','Best')
ylabel('F-score')
set(gca,'FontSize', 12);
%% Figure 7, overcluster
test_k = [3 4 5 6 7 8 9 10 15 20 25 30 35 40 45 50 75 100 150];
outtable = cell(size(test_k));
nmi = zeros(size(test_k));
m = cell(size(test_k));
parfor i = 1:length(test_k)
    rng(42)
    km_l = kmeans(data,test_k(i));%,'Options',opts);
    m{i} = C2G(data,km_l,km_l);
    [outtable{i},nmi(i)] = m{i}.show_f_score(km_l);
end
fscores = cellfun(@(a) a(:,6),outtable,'UniformOutput',false)';
save overcluster.mat
% groups = cellfun(@(a) ones(size(a)),fscores,'UniformOutput',false);
% for i = 1:length(test_k)
%     groups{i} = groups{i} * test_k(i);
% end
figure
plot(nmi,'-^','Color', [0 0 0])
hold on;
plot(cellfun(@nanmean,fscores),'-o','Color', [0.5 0.5 0.5])
axis([0 length(test_k)+1 0 1]);
legend('NMI','Mean F-score')
xlabel('K');
ylabel('Scores');
set(gca,'XTick',1:length(test_k));
set(gca,'XTickLabel',test_k);
set(gca,'FontSize',12);
%% Figure 8, outliers
noise = [0 1 2 4 8 16 32];
exclue_noise = [0 1 2 4 10 20 50]./100;
outtable = cell(size(noise));
nmi = zeros(size(noise));
for i = 1:length(noise)
    rng(42)
    noise_l = make_noise(ori_l, noise(i));%,'Options',opts);
    label = cluster_ungated(data,noise_l,'outliers', exclue_noise(i));
    m{i} = C2G(data,label,noise_l,'outliers',exclue_noise(i));
    [outtable{i},nmi(i)] = m{i}.show_f_score(ori_l);
end
save outliers.mat
figure;
plot(1:length(noise), nmi, '-^','Color', [0 0 0]);
hold on 
f = cellfun(@(x) mean(x(:,6)),outtable);
plot(1:length(noise), f, '-o','Color', [0.5 0.5 0.5]);
axis([0.5 length(noise)+1 0 1])
set(gca,'XTick',1:length(noise));
set(gca,'XTickLabel',noise);
xlabel('Outlier Percentage')
ylabel('Score')
legend({'NMI','Mean F-score'},'Location','best')
set(gca, 'FontSize', 12);
%% Figure S2, Table S2, phenograph 
load('testdata/Phenograph_labels.mat','labels');
m_phno = C2G(data, labels, labels, 'markernames',markers,'trivial_gate',200);
m_phno.show_f_score(labels);
m_phno.view_gates_contour(data,markers,'ignore_small',0);
% phno_flowtype = compareflowtype(data, labels);

%% Figure S3, Table S3, tSNE
load testdata/tSNE_label.mat;
m_tsne = C2G(data, tSNE_label, tSNE_label, 'markernames',markers, 'trivial_gate',200);
m_tsne.show_f_score(tSNE_label);
m_tsne.view_gates(data,markers,'ignore_small',0, 'n_line',2,'onepanel',true);

%% Figure S4, Table S4, flowSOM
d = csvread('build/data.csv')';
net = selforgmap([3 3]);
net = train(net, d);
[~, l] = max(net(d));
l = l';
d = d';

m_som = C2G(d, l, l, 'markernames', markers);
out = m_som.show_f_score(l);
m_som.view_gates(d, markers);
m_som.view_gates_contour(d, markers);
%% Figure S5, Table S5
d = csvread('build/data.csv');
l = csvread('build/label_flowmeans.csv',0,0);

m_fm = C2G(d, l, l, 'markernames', markers);
out = m_fm.show_f_score(l);
m_fm.view_gates(d, markers);
m_fm.view_gates_contour(d, markers);
%% Compare to Random shuffle
rng(9464)
km_l = kmeans(data,10);
new_km_l = km_l; % Since all populiations is known, no need to pre-cluster
% Perform the anlysis
m_km = C2G(data,km_l,km_l,'markernames',markers,'trivial_gate',200,'ratio_trivial_gate',0.1);
% Visualize the results
m_km.view_gates(data,markers,'n_lines',4,'ignore_small',300);
m_km.show_f_score(km_l);

maxdepth = max(cell2mat(m_km.depth));

rng(42)
nmi = zeros(500,1);
for i = 1:500
    m_rnd = C2G(data,km_l,km_l,'markernames',markers,'trivial_gate',200,...
        'ratio_trivial_gate',0.1,'maxdepth',maxdepth,'randpair',true);
    [~,nmi(i)] = m_rnd.show_f_score(km_l);
end
figure
histogram(nmi);
hold on;
plot([0.77 0.77],[0 80],'r')
xlabel('NMI')
ylabel('Frequency')
set(gca,'FontSize',12);
%% Larger data
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
m.view_gates_contour(d,marker(marker_idx),'ignore_small',3000,'n_lines',4);
warning(w);
m.show_f_score(label);
%% Timing C2G
% Impact by number of target populations
clear
[ori_data, marker] = readfcs_v2('testdata/bigdata/TIN_BLD1_Untreated_Day3.fcs');
data = flow_arcsinh(ori_data,5);
marker_idx = [10,16,17,18,19,21,22,24,29,31,32,33,39,40,41,46,47,49,50,51,52];
d = data(marker_idx,:)';
markers = marker(marker_idx);
rng(9464)
label = kmeans(d,10);
run_time = cell(91,1);
parfor k = 0:90
    i = floor(k/10) + 1 ;
    j = mod(k, 10) + 1;
    m_idx = j:(j+i-1);
    m_idx(m_idx>10) = m_idx(m_idx>10) - 10;
    ln = label;
    ln(~ismember(label, m_idx)) = 0;
    cln = cluster_ungated(d,ln);
    tic;C2G(d, cln, ln,'markernames',markers,'trivial_gate',200,'showdetail',0);t = toc;
    run_time{k+1} = t;
end
run_time = cell2mat(run_time);
save timeall.mat`
pop = repelem(1:10, 10);
boxplot(run_time, pop(1:91));
axis([0.5 10.5 0 9000])
xlabel('Number of Target Populations')
ylabel('Running Time/s')
% Impact by number of cells in the first iteration
cellnum = [20000 20000 40000 40000 60000 60000 80000 80000 100000 100000 120000 120000 140000];
run_time_cellnum = cell(size(cellnum));
rng(42);
parfor i = 1:length(cellnum)
    l = label;
    cidx  = randsample(size(d,1), cellnum(i));
    subd = d(cidx,:);
    subl = l(cidx);
    tic;C2G(subd, subl, subl,'trivial_gate',200,'maxdepth',1,'showdetail',0);t = toc;
    run_time_cellnum{i} = t;
end
run_time_cellnum = cell2mat(run_time_cellnum);
% Impact by number of protein markers in the first iteration
dimnum = [3 6 9 12 15 18 21];
run_time_dimnum = cell(31,1);
rng(42);
parfor k = 1:31
    i = floor(k/5) + 1 ;
    j = mod(k, 5) + 1;
    dn = dimnum;
    m_idx = j:(j+dn(i)-1);
    m_idx(m_idx>21) = m_idx(m_idx>21) - 21;
    subd = d;
    subd = subd(:,m_idx);
    tic;C2G(subd, label, label,'trivial_gate',200,'maxdepth',1,'showdetail',0);t = toc;
    run_time_dimnum{k} = t;
end
run_time_dimnum = cell2mat(run_time_dimnum);
dimnums = repelem(dimnum, 5);
save time.mat
fontsize = 12;
figure('pos',[100 100 800 300]);
subplot(1,2,1)
boxplot(run_time_cellnum, cellnum);
% axis([0 145000 0 90]);
xlabel('Number of Cells')
ylabel('First iteration Time/s')
set(gca,'FontSize',fontsize);
xtickangle(30)
subplot(1,2,2)
boxplot(sqrt(run_time_dimnum), dimnums(1:31));
% axis([1 21.5 0 9]);
xlabel('Number of Protein Markers')
ylabel('$\bf{First\enspace iteration\enspace \sqrt{Time}/ \sqrt{s}}$','Interpreter','latex');
set(gca,'FontSize',fontsize-2);



% Impact by number of known populations
run_time_knowpop = cell(91,1);    
rng(42);
parfor k = 1:91
    i = floor(k/10) + 1 ;
    j = mod(k, 10) + 1;
    m_idx = j:(j+i-1);
    m_idx(m_idx>10) = m_idx(m_idx>10) - 10;
    subl = label;
    subl(~ismember(label, m_idx)) = 0;
    %subl = kmeans(d, knownpop(i));
    subd = d;
    tic;C2G(subd, subl, subl,'trivial_gate',200,'maxdepth',1);t = toc;
    run_time_knownpop{k} = t;
end
save timeall_1pop.mat
run_time_knowpop = cell2mat(run_time_knownpop);
knownpop = repelem(1:10, 10);
boxplot(run_time_knowpop, knownpop(1:91));
figure('pos',[100 100 800 600]);
fontsize = 12;
subplot(2,2,1)
scatter(knownmarker, run_time, 15, [0  0 0],'filled');
axis([0.5 10.5 0 5000]);
xlabel('Number of Target Populations')
ylabel('Time/s')
set(gca,'FontSize',fontsize);
subplot(2,2,2)
scatter(knownpop, run_time_knownpop, 15, [1 0 0],'filled');
axis([0.5 10.5 0 110]);
xlabel('Number of Target Populations')
ylabel('First iteration Time/s')
set(gca,'FontSize',fontsize);

figure('pos',[100 100 800 300]);
subplot(1,2,1)
boxplot(run_time_cellnum, cellnum);
% axis([0 145000 0 90]);
xlabel('Number of Cells')
ylabel('First iteration Time/s')
set(gca,'FontSize',fontsize);
xtickangle(30)
subplot(1,2,2)
boxplot(sqrt(run_time_dimnum), dimnum);
% axis([1 21.5 0 9]);
xlabel('Number of Protein Markers')
ylabel('$\bf{First\enspace iteration\enspace \sqrt{Time}/ \sqrt{s}}$','Interpreter','latex');
set(gca,'FontSize',fontsize-2);
%% Pool multiple samples
[ori_d,~,markers] = load_mul_fcs('testdata/bigdata');
% Select protein markers
marker_idx = [10,16,17,18,19,21,22,24,29,31,32,33,39,40,41,46,47,49,50,51,52];
d = ori_d(:,marker_idx);
rng(42);
label = kmeans(d,10);
tic;m = C2G(d, label, label,'markernames',markers(marker_idx),'trivial_gate',2000,'ratio_trivial_gate',0.1);toc;
m.show_f_score(label);
w = warning ('off','all');
m.view_gates(d,markers(marker_idx),'ignore_small',3000,'n_lines',4);
warning(w);
% Elapsed time is 2562.297240 seconds

%% Compare to heatmap
rng(1);

figure;
subplot(1,3,1)
a = mvnrnd([0 2], [0.5 0;0 5]./7, 500);
b = mvnrnd([1.5 1.5], [1 0.8;0.8 1]./3, 2500);
c = mvnrnd([2 0], [5 0;0 0.5]./7, 1500);
hold on
scatter(a(:,1), a(:,2),5, 'g','.');
scatter(c(:,1), c(:,2),5, 'b','.');
scatter(b(:,1), b(:,2),5, 'r','.');
axis([-2 5 -2 5])
xlabel('Marker 1')
ylabel('Marker 2')
subplot(1,3,2)
a = mvnrnd([0 2], [1 0;0 1]./7, 500);
b = mvnrnd([1.5 1.5], [1 0;0 1]./7, 2500);
c = mvnrnd([2 0], [5 0;0 0.5]./3, 1500);
hold on
scatter(a(:,1), a(:,2),5, 'g','.');
scatter(c(:,1), c(:,2),5, 'b','.');
scatter(b(:,1), b(:,2),5, 'r','.');
axis([-2 5 -2 5])
xlabel('Marker 1')
ylabel('Marker 2')
subplot(1,3,3)
h = heatmap({'Green', 'Blue', 'Red'},{'Marker 1', 'Marker 2'},[0 2; 1.5 1.5; 2 0]' );
h.Colormap = jet;


%% Show why need overlapping gates
rng(1);
a = mvnrnd([0 0], eye(2)/10, 1000);
b = mvnrnd([0.2 0], eye(2)/10, 1000);
d = [a;b];
c = [repmat([1,0,0],size(a,1),1);repmat([0,0,1],size(b,1),1)];
figure;
subplot(1,2,1)
scatter(d(:,1),d(:,2),10,c);
hold on
plot(b(convhull(b),1),b(convhull(b),2), '-b',  'LineWidth', 3)
plot(a(convhull(a),1),a(convhull(a),2), '-r',  'LineWidth', 3)
title('Know One Gate')
axis([-1 1.5 -1 1.5])
xlabel('Marker 1')
ylabel('Marker 2')
subplot(1,2,2)
scatter(d(:,1),d(:,2),10,c);
hold on
plot(d(convhull(d),1),d(convhull(d),2), '-k', 'LineWidth', 3)
title('Know All Gate')
axis([-1 1.5 -1 1.5])
xlabel('Marker 1')
ylabel('Marker 2')


figure;
subplot(1,4,1)
scatter3(d(:,1), d(:,2), d(:,3), 10, c, 'Marker','.')
subplot(1,4,2)
scatter(d(:,1),d(:,2), 10,c, 'Marker','.')
subplot(1,4,3)
scatter(d(:,2),d(:,3), 10,c, 'Marker','.')
subplot(1,4,4)
scatter(d(:,1),d(:,3), 10,c, 'Marker','.')


%figure('pos',[100 100 400 400])
[x,y,z] = sphere;
subplot(1,4,1)
hold on 
surf(x,y,z, ones(21,21)) % centered at (3,-2,0) 
surf(x+1.2,y+1.2,z+1.2, ones(21,21)*0.5) % centered at (0,1,-3)
xlabel('Dim 1')
ylabel('Dim 2')
zlabel('Dim 3')
subplot(1,4,2)
hold on 
surf(x,y,z, ones(21,21)) % centered at (3,-2,0) 
surf(x+1.2,y+1.2,z+1.2, ones(21,21)*0.5) % centered at (0,1,-3)
xlabel('Dim 1')
ylabel('Dim 2')
zlabel('Dim 3')

subplot(1,4,3)
hold on 
surf(x,y,z, ones(21,21)) % centered at (3,-2,0) 
surf(x+1.2,y+1.2,z+1.2, ones(21,21)*0.5) % centered at (0,1,-3)
xlabel('Dim 1')
ylabel('Dim 2')
zlabel('Dim 3')

subplot(1,4,4)
hold on 
surf(x,y,z, ones(21,21)) % centered at (3,-2,0) 
surf(x+1.2,y+1.2,z+1.2, ones(21,21)*0.5) % centered at (0,1,-3)
xlabel('Dim 1')
ylabel('Dim 2')
zlabel('Dim 3')

%% Rare cell population
[~, marker, ~,~,ori_d] = readfcs_v2('testdata/rare population/108.fcs');
d = flow_arcsinh(ori_d, 100)';
l = csvread('testdata/rare population/108.csv');

rng(9464)
label = cluster_ungated(d,l);
% Perform the anlysis
m = C2G(d,label,l,'markernames',marker);
% Visualize the results
m.view_gates(d,marker,'n_lines',1,'ignore_small',0,'onepanel',true);
% m.view_gates_contour(d,marker,'n_lines',3,'ignore_small',0,'onepanel',true);
m.show_f_score(l);


[~, ori_marker, ~,~,ori_d] = readfcs_v2('testdata/rare population/visne_mrd.fcs');
fd = flow_arcsinh(ori_d, 5)';
subidx = [9 12 13 15 16 17 18 19 21 24:30 32 34 36 37 38 39 42];
d = fd(:,subidx);
marker = ori_marker{subidx};

net = selforgmap([10 10]);
net = train(net, d');
[~, l] = max(net(d'));
l = l';



[ori_d, ori_marker, ~,~] = readfcs_v2('testdata/rare population/Mosmann_rare_notransform.fcs');
fd = flow_arcsinh(ori_d, 100)';
l = csvread('testdata/rare population/X-shift-label.csv');
l = l-15;

subidx = [7 8 9 11 12 13 14 15 16 17 18 19 20];
d = fd(:,subidx);
marker = ori_marker(subidx);

l(l~=4) = 0;
l(l==4) = 1;
label = cluster_ungated(d,l);
m = C2G(d,label,l,'markernames',marker);
% Visualize the results
m.view_gates(d,marker,'n_lines',1,'ignore_small',0,'onepanel',true);
% m.view_gates_contour(d,marker,'n_lines',3,'ignore_small',0,'onepanel',true);
m.show_f_score(l);
