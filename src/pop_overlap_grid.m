function adj = pop_overlap_grid(numc)
np = size(numc,2)-1;% Ignore first column (outliers)
adj = eye(np);
for i = 1:np-1
    for j = i+1:np
        adj(j,i) = sum(sqrt(numc(:,i+1)/sum(numc(:,i+1)).*numc(:,j+1)/sum(numc(:,j+1))));
%         overlap_idx = numc(:,i+1).*numc(:,j+1)>0;
%         adj(j,i) = max(sum(numc(overlap_idx,i+1))/sum(numc(:,i+1)),sum(numc(overlap_idx,j+1))/sum(numc(:,j+1)));
        adj(i,j) = adj(j,i);
    end
end
%areas = areas + areas' - eye(length(all_labels));
adj(isnan(adj))=0;
%adj(isnan(areas))=0;

