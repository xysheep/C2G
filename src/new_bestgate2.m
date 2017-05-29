function [penalty, group, children, flag_end] = new_bestgate2(label,ori_label,means,covs,d12,overlap2d)


%% Overlap between different populations
% Implement this part in C++/ vectorization 
%[adj,~,all_labels] = pop_overlap(x,y,label,density,ig_ratio); 
[adj,overlap] = gaussian_pop_overlap(means,covs,label,d12); 
overlap = overlap + overlap2d(label,label,d12(1),d12(2))*100;
adj = adj.*(~overlap2d(label,label,d12(1),d12(2)));
if length(label)>2
    group = cluster(overlap);
elseif length(label)==2
    group = [1;2];
elseif length(label)==1
    group = 1;
end
children = cell(length(label),1);
penalty = sum(sum(triu(adj)));
for i = 1:max(group)
    if sum(ismember(label(group==i), ori_label))>0
        children{i} = label(group==i);
    end
    penalty = penalty - sum(sum(triu(adj(group==i))));
end
children = children(~cellfun('isempty',children));

flag_end = max(group)==1;

