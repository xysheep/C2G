function visulize_gating_hierechy(d,M,l,ori_label,markers)
n_markers = size(d,2);
% gated_d = d(l~=0,:);
% gated_l = l(l~=0);
gated_d = d;
gated_l = l;
%num_c = length(unique(gated_l))+1;
if ~exist('markers','var')
    markers = strread(num2str(1:n_markers),'%s');
end
num_c = length(unique(gated_l(gated_l~=0)))+1;
k = 1;
pairorder = zeros(n_markers);
for i=1:n_markers-1
    for j=i+1:n_markers
        k = k+1;
        pairorder(i,j) = k; 
    end
end
n_graph = nnz(M);
[i,j]=find(M);

% This section decide the order of 2D graph shown
o=zeros(size(i));
queue = CQueue();
queue.push(find(i==1,1));
k=1;
while k<=length(i)
    o(k) = queue.pop();
    for ele = find(i==j(o(k)))'
        queue.push(ele)
    end
    k = k + 1;
end
i=i(o);
j=j(o);

% Draw the figure
% To make the color consistent, pre-assign different color to each population 
colorset = varycolor(num_c);
if exist('ori_label','var')
    unused_idx = unique(l(ori_label==0));
    t = 1:num_c;
    used_idx = t(~ismember(t,unused_idx));
    colorset(unused_idx,:) = repmat([0.7 0.7 0.7],length(unused_idx),1);
    colorset(used_idx,:) = varycolor(num_c-length(unused_idx));
end
num2draw = 100000;
for k=1:n_graph
    subplot(ceil(sqrt(n_graph)),ceil(sqrt(n_graph)),k)
    clusters2draw = decode(M(i(k),j(k)),num_c);
    % This visulzation step didn't consider which cells are carried in this
    % step. Need to rewrite a lot here.
    [dim1,dim2] = find(pairorder==j(k));
    idx = find(ismember(gated_l,clusters2draw));
    if length(idx)>num2draw
        rng(9464);
        y = randsample(length(idx),num2draw);
        idx = idx(y);
    end
    visulizeGroups_colorset(gated_d(idx,dim1),gated_d(idx,dim2),gated_l(idx),colorset);
    title(sprintf('%d',length(unique(clusters2draw))));
    xlabel(markers{dim1})
    ylabel(markers{dim2})
    if k==13
        disp(unique(clusters2draw));
    end
end
    



