function [downgroup,A]=cluster(adj)
A=adj;
downgroup = ones(size(A,1),1);
queue=unique(downgroup);
while size(queue)>0
    g=find(downgroup==queue(1));
    B=modularity(A,g);
    [eigen1,~]=eigs(B,1);
    new_group_idx=(eigen1>0);
    if sum(new_group_idx)>1 && sum(new_group_idx)<length(g)
        downgroup(g(new_group_idx))=max(downgroup)+1;
        queue=[queue max(downgroup)];
    else
        queue=queue(2:end);
    end
end
