function main(datafile, clusterfile, varargin)
d = csvread(datafile);
l = csvread(clusterfile);
if size(d,1) ~= length(l)
    fprintf("Data and labels have differnt number of cells\n");
    return;
end
numoption = {'ratio_trivial_gate', 'trivial_gate', 'grid_size','randpair','maxdepth','outliers','fontsize','n_lines','ignore_small','onepanel'};
c2g_option = { 'ratio_trivial_gate', 'trivial_gate','markernames','color','showdetail', 'grid_size','randpair','maxdepth','outliers'};
viewgates_option = {'fontsize','n_lines','ignore_small','onepanel'};
c2gvar = cell(0);
viewvar = cell(0);
for i = 1:2:size(varargin,2)
    if ismember(varargin{i}, numoption) && ischar(varargin{i+1})
        varargin{i+1} = str2double(varargin{i+1});
    end
    if ismember(varargin{i}, c2g_option)
        c2gvar{end+1} = varargin{i};
        c2gvar{end+1} = varargin{i+1};
    end
    if ismember(varargin{i}, viewgates_option)
        viewvar{end+1} = varargin{i};
        viewvar{end+1} = varargin{i+1};
    end
end
label = cluster_ungated(d,l);
[m,~,markernames] = C2G(d, label, l, c2gvar{:});
m.view_gates(d, markernames,  viewvar{:});
m.show_f_score(l);