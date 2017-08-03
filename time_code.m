load testdata\Screen_PACA-1-3_008.mat
ori_data = data;
data = ori_data(7:16,ori_data(21,:)>0.1)';
markers = marker_names(7:16);

rng(9464)
km_l = kmeans(data,10);
km_l(km_l > 6) = 0;
% Precluster the ungated cells
label = cluster_ungated(data,km_l);
new_km_l = km_l;
% Precompute the local density
[~, ~, density] = compute_density(data,label);
% Perform the anlysis
m = C2G(data,label,new_km_l,density,'ignore_ratio',0.2);
m.view_gates(data,markers,20);
m.show_f_score(new_km_l);