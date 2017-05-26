rng(9464)
km_l = kmeans(data,10);
% Precluster the ungated cells
label = cluster_ungated_gmm(data,km_l);
new_km_l = km_l;
% Precompute the local density
[means, covs, density] = compute_density(data,label);
% Perform the anlysis
m = C2G(data,label,new_km_l,means,covs,density,'ignore_ratio',0.2);
m.visulize_gating_sequence(data,markers,4,100);
m.show_f_score(new_km_l);