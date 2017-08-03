function [gated_labels,main_members,ignore_perc,boundaries,flag,adj] = new_bestgate(x,y,label,ori_label,ori_main_members,density,ig_ratio)


%% Overlap between different populations
% Implement this part in C++/ vectorization 
[adj,all_labels] = pop_overlap(x,y,label,density,ig_ratio,ori_label); 
%[adj,overlap,all_labels] = gaussian_pop_overlap(means,covs,label,d12); 

%adj;
m = mcl(adj);
% figure;imagesc(adj);
% figure;imagesc(m);
%sum(sum(m,2)>=1)
m = round(m,3);
% m = cluster(-log(adj+max(0.0001,min(adj(:)))*0.01));
gates = zeros(size(label));
for i = 1:length(all_labels)
    [~,gates(label==all_labels(i))] = max(m(:,i));
%     gates(label==all_labels(i)) = m(i);
    % cell and gate label are not necessary related. They are just index of
    % different things
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
uniq_gate = unique(gates);
uniq_gate(uniq_gate==0) = [];
%fprintf('clusters = %2d\t',length(uniq_gate));
%tabulate(m(:));
if length(uniq_gate)>1
    gated_labels = cell(1,length(uniq_gate));
    main_members = cell(1,length(uniq_gate));
    ignore_perc = cell(1,length(uniq_gate));
    boundaries = cell(1,length(uniq_gate));
    for i = 1:length(uniq_gate)
        uniq_label = unique(label(gates==uniq_gate(i)));
        uniq_ori_label = uniq_label(ismember(uniq_label,ori_label));
        if ~isempty(uniq_ori_label)   % gates only from ungated cells should be ignored
            %label_in_one_polygon = cell(1,length(uniq_ori_label));
            label_in_one_polygon = cell(1,length(uniq_label));
            for j = 1:length(uniq_label)
                %idx_in_2D = (gates == uniq_gate(i) & label == uniq_ori_label(j));
                idx_in_2D = (gates == uniq_gate(i) & label == uniq_label(j));
                label_in_one_polygon{j} = idx_in_2D;
            end
            [gated_labels{i},ignore_perc{i},boundaries{i}] = gated_cells(x,y,label_in_one_polygon,ismember(label,uniq_label),density);
            main_members{i} = uniq_label(ismember(uniq_label,ori_main_members));
        end
    end
    main_members = main_members(~cellfun('isempty',gated_labels));
    gated_labels = gated_labels(~cellfun('isempty',gated_labels));
    ignore_perc = ignore_perc(~cellfun('isempty',ignore_perc));
    boundaries = boundaries(~cellfun('isempty',boundaries));
else
    main_members = {ori_main_members};
    gated_labels = {1:length(x)};
    ignore_perc = 0;
    [bx,by] = findboundary(x,y,x,y);
    boundaries = {[bx by]};
    flag = false;
end




