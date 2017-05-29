function T = C2G2(d,l,ori_l,means,covs,overlap2d,varargin)%ig_ratio,markernames,col)
% C2G perform the analysis return a gatingTree object that store the
% obtained gating hierarchy.
%       T = C2G(d,l,ori_l,...) "d" is the M-by-N data matrix where M is the
%       number of markers. "l" and "ori_l" are M-by-1 matrix represent cell
%       labels after and before pre-cluster. 
%       T = C2G(d,l,ori_l,density,...) "density" is the precomputed local
%       density. It can be computed outside C2G use function
%       "compute_density". If this parameter is not provided, C2G will
%       compute local density later.
% Optional parameter:
%       ignore_ratio: percentage of low density cells ignored when compute
%       overlap between different populations. If the data is expected to
%       be highly noisy, increase this value for better performance.
%       Default is 0.05. 

n_markers = size(d,2);
% Compute local density

% Initiate other parameters
pnames = { 'ignore_ratio','markernames','color'};
dflts  = { 0.05          ,[]           ,[]};
[~,markernames,col] = internal.stats.parseArgs(pnames,dflts,varargin{:});


queue = CQueue();
queue.push(1);
l = unique(l);
ori_l = unique(ori_l);
T = gatingTree(l,4);
step = 0;
% T is the object representing the tree of gating hierachy
while ~queue.isempty()
    step = step + 1;
    node_id = queue.pop();
        best_penalty = inf;
        best_children = cell(0);
        best_pair = 0;
        sub_l = T.cell_label{node_id};
        sub_ori_l = ori_l(ismember(ori_l,sub_l));
        for i = 1:n_markers - 1
            for j = i+1:n_markers
                
                d12 = [i j];
                [penalty, ~, children, flag_end] = new_bestgate2(sub_l,unique(sub_ori_l),means,covs,d12,overlap2d);
                              
                if penalty < best_penalty && ~flag_end
                    best_penalty = penalty;
                    best_children = children;
                    best_pair = d12;
                end
            end
        end

        
        if ~isempty(best_children)
            T.setdim(node_id,best_pair);
            for i_gate=1:length(best_children)
                T.addnode(node_id, best_children{i_gate});
                if length(best_children{i_gate})>1
                    queue.push(T.numNode);
                end
            end
        end
end


    