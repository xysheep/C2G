function [gated_labels,main_members,boundaries,flag,adj] = new_bestgate_grid(x,y,label,ori_label,ori_main_members,GRID_SIZE)
all_labels = unique(label);
all_labels(all_labels==0 | histc(label,all_labels)<1) = [];
ori_label(ori_label==0) = [];
%% Build Grid
x_inv = range(x)/GRID_SIZE+0.0001;
y_inv = range(y)/GRID_SIZE+0.0001;
nx = floor((x-min(x))/x_inv);
ny = floor((y-min(y))/y_inv);
code = nx*100 + ny;
[grid,~,ic] = unique(code);
% numc = zeros(size(grid,1),1+length(all_labels));
% for i = 1:length(nx)
%     if label(i) == 0
%         numc(grid==code(i),1) = numc(grid==code(i),1) + 1;
%     else
%         numc(grid==code(i),find(all_labels==label(i))+1) = numc(grid==code(i),find(all_labels==label(i))+1) + 1;% Ungated 
%     end
%     % This step is to compute number of each cell of each cell population in each grid 
% end
numc = zeros(size(grid,1),1+length(all_labels));
counter  =  1;
for t = [0,all_labels']
    lt = label == t;
    if sum(lt) > 0
        numc(:,counter) = histc(ic(lt), 1:length(grid));
    end
    counter = counter + 1;
end
%% Overlap between different populations
% Implement this part in C++/ vectorization 
%[adj,all_labels] = pop_overlap(x,y,label,density,ig_ratio); 
adj = pop_overlap_grid(numc); % The first column is ignored
% figure;imagesc(adj1),figure;imagesc(adj);
% figure;imagesc(mcl(adj1)),figure;imagesc(mcl(adj));
%[adj,overlap,all_labels] = gaussian_pop_overlap(means,covs,label,d12); 

%adj;
m = mcl(adj);
% figure;imagesc(adj);
% figure;imagesc(m);
%sum(sum(m,2)>=1)
m = round(m.*1000)./1000;
% m = cluster(-log(adj+max(0.0001,min(adj(:)))*0.01));
% gates = zeros(size(label));
% for i = 1:length(all_labels)
%     [~,gates(label==all_labels(i))] = max(m(:,i));
% %     gates(label==all_labels(i)) = m(i);
%     % cell and gate label are not necessary related. They are just index of
%     % different things
% end
g = zeros(size(all_labels));
for i = 1:length(all_labels)
    [~,g(i)] = max(m(:,i));
end
%% If only one cluster, remove the cluster with smallest clustering coefficient and do mcl again

% % betweenness or closeness may be better metric to remove noise ~ Do it tomorrow 
% if length(unique(gates))==1
%     cc = clustering_coefficients(adj);
%     [~,rm_idx] = min(cc);
%     rm_label = all_labels(rm_idx);
%     all_labels(rm_idx) = [];
%     adj(rm_idx,:) = [];
%     adj(:,rm_idx) = [];
%     m = mcl(adj);
%     gates = zeros(size(label));
%     for i = 1:length(all_labels)
%         [~,gates(label==all_labels(i))] = max(m(:,i));
%         % cell and gate label are not necessary related. They are just index of
%         % different things
%     end
%     gates(label == rm_label) = length(cc);
% end

%% Find all cells in each gate
flag = true;
uniq_gate = unique(g);
%uniq_gate(uniq_gate==0) = [];
%fprintf('clusters = %2d\t',length(uniq_gate));
%tabulate(m(:));
if length(uniq_gate)>1
    gated_labels = cell(1,length(uniq_gate));
    main_members = cell(1,length(uniq_gate));
    boundaries = cell(1,length(uniq_gate));
    for i = 1:length(uniq_gate)
        uniq_label = find(g==uniq_gate(i));
        uniq_ori_label = uniq_label(ismember(all_labels(uniq_label),ori_label));
        if ~isempty(uniq_ori_label)   % gates only from ungated cells should be ignored
            [gated_labels{i},grid_boundary] = gated_cells_grid(grid,numc,ic,uniq_label);
            main_members{i} = all_labels(uniq_label(ismember(all_labels(uniq_label),ori_main_members)));
            if size(grid_boundary,2)==2
                grid_x = grid_boundary(:,2);
                grid_y = grid_boundary(:,1);
                rm_idx = isnan(grid_x)|isnan(grid_y);
                grid_x(rm_idx) = [];
                grid_y(rm_idx) = [];
                tmpx = [grid_x;grid_x;grid_x+1;grid_x+1];
                tmpy = [grid_y;grid_y+1;grid_y;grid_y+1];
                %[bx,by] = findboundary_grid(x,y);
                [bx,by] = findboundary_grid(tmpx,tmpy);
%                 boundaries{i} = [min(x) + x_inv*grid_x min(y) + y_inv*grid_y ];
                boundaries{i} = [min(x) + x_inv*bx min(y) + y_inv*by ];
            end
        end
    end
    
    main_members = main_members(~cellfun('isempty',boundaries));
    gated_labels = gated_labels(~cellfun('isempty',boundaries));
    boundaries = boundaries(~cellfun('isempty',boundaries));
else
    main_members = {ori_main_members};
    gated_labels = {1:length(x)};
    [bx,by] = findboundary(x,y,x,y);
    boundaries = {[bx by]};
    flag = false;
end




