function main(fcsfile, clusterfile, cofactor)
[~,markernames, ~,~, ori_d] =readfcs_v2(fcsfile);
d = flow_arcsinh(ori_d,cofactor)';
l =csvread(clusterfile);
if size(d,1) ~= length(l)
    fprintf("FCS and Labels have differnt number of cells\n");
    return;
end


label = cluster_ungated(d,l);
m = C2G(d, label, l, 'markernames', markernames);
m.view_gates(d, markernames, 'ignore_small', 0, 'n_lines',4);
m.show_f_score(l);