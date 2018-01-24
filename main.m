function main(datafile, clusterfile, varargin)
%% C2G is developed to generate nested 2D gating hierarchy from clustered flow and mass cytometry data.
% Required Arguments:
%       datafile: Path to a CSV file of data matrix. This file should have a
%               n-by-m matrix without header. n is the number of cells and
%               m is number of markers. This data should already be
%               compensated and transformed.
%       clusterfile: Path a CSV file of cluster labels. This file
%               should have a n-by-1 matrix without header. n is the number
%               of cells. Each row correspondong to one cell in datafile. 0
%               means unlabeled.
% Optional Arguments:
%       Details of optional arguments are avaliable at https://github.com/xysheep/C2G
if ~exist('datafile','var') || ~exist('clusterfile','var')
    fprintf(['C2G is developed to generate nested 2D gating hierarchy from clustered flow and mass cytometry data.\n',...
'    Required Arguments:\n',...
'          datafile: Path to a CSV file of data matrix. This file should have a\n',...
'                  n-by-m matrix without header. n is the number of cells and\n',...
'                  m is number of markers. This data should already be\n',...
'                  compensated and transformed.\n',...
'          clusterfile: Path a CSV file of cluster labels. This file\n',...
'                  should have a n-by-1 matrix without header. n is the number\n',...
'                  of cells. Each row correspondong to one cell in datafile. 0\n',...
'                  means unlabeled.\n',...
'    Optional Arguments:\n',...
'          Details of optional arguments are avaliable at https://github.com/xysheep/C2G\n'])
    return
end


d = csvread(datafile);  
l = csvread(clusterfile);
if size(d,1) ~= length(l)
    fprintf('Data and labels have differnt number of cells\n');
    return;
end
numoption = {'maximum_cluster', 'ignore_small_bin','cluster_amount','outliers','visualize','ratio_trivial_gate', 'trivial_gate', 'grid_size','randpair','maxdepth','outliers','fontsize','n_lines','ignore_small','onepanel'};
c2g_option = { 'ratio_trivial_gate', 'trivial_gate','markernames','color','showdetail', 'grid_size','randpair','maxdepth','outliers'};
viewgates_option = {'fontsize','n_lines','ignore_small','onepanel'};
precluster_option = {'maximum_cluster', 'ignore_small_bin','cluster_amount','outliers','visualize'};
c2gvar = cell(0);
viewvar = cell(0);
precvar = cell(0);
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
    if ismember(varargin{i}, precluster_option)
        precvar{end+1} = varargin{i};
        precvar{end+1} = varargin{i+1};
    end
end

pnames = { 'cofactor'};
dflts  = { 0};
[cofact] = internal.stats.parseArgs(pnames,dflts,varargin{:});
if ischar(cofact)
    cofact = str2double(cofact);
end
if cofact > 0
    d = flow_arcsinh(d', cofact)';
end
label = cluster_ungated(d,l,precvar{:});
[m,~,markernames] = C2G(d, label, l, c2gvar{:});
m.view_gates(d, markernames,  viewvar{:});
m.view_gates_contour(d, markernames,  viewvar{:});
m.show_f_score(l);