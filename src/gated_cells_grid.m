function [best_idx,best_boundary] = gated_cells_grid(grid,numc,ic,target_labels)
% This function find all cells in the union of polygons of all gated
% populations. 

%% Testing command:
%%
% Use the following scripts to test this file
%   
%   fdname = '../FlowData/Primary Murine T cell Data';
%   [data,labels]=load_mul_fcs(fdname,'ctr.fcs');
%	i=8;j=27;
%   gatelabels = {find(labels==1) find(labels==2)};
%   gated_cells(data(:,i),data(:,j),gatelabels);


%% Main part
max_fscore = 0;
best_idx = [];
best_boundary = 0;
x = mod(grid,100);
y = floor(grid/100);
percs = [0.0 0.05 0.1 0.2 0.3];
for ignore_perc = percs
    idx_cell = cell(length(target_labels),1);
    boundariesx = [];
    boundariesy = [];
    for i = 1:length(target_labels)
        n_cell = sum(numc(:,1+target_labels(i)));
        if ( n_cell > 20) && ignore_perc >= 0.05
            [s,si] = sort(numc(:,1+target_labels(i)));
            idx = find(cumsum(s)>=sum(s)*ignore_perc + 1,1);
            idx = si(idx:end); 
            if isempty(idx)
                idx = si(end);
            end
        else
            [s,si] = sort(numc(:,1+target_labels(i)));
            idx = find(cumsum(s)>=0 + 1,1);
            idx = si(idx:end); 
            if isempty(idx)
                idx = si(end);
            end
        end
        xg = x(idx);
        yg = y(idx);
        [boundaryx,boundaryy] = findboundary_grid(xg,yg);
        if length(boundaryx) > 2 % At least 3 points to be a polygon
            idx_cell{i} = find(inpolygon_grid(x,y,boundaryx,boundaryy));
            if isempty(boundariesy)
                boundariesx = boundaryx;
                boundariesy = boundaryy;
            else
                try
                    [boundariesx,boundariesy] = polybool('union',boundariesx,boundariesy,boundaryx,boundaryy);
                catch
                    1;
                end
            end
        end
    end
    gated_idx = unique(cell2mat(idx_cell));

    f_score = fscore_grid(numc,target_labels,gated_idx);
    if f_score > max_fscore
        max_fscore = f_score;
        best_idx = gated_idx;
        best_boundary = [boundariesx boundariesy];
    end
end
best_idx = find(ismember(ic,best_idx));
%% Visulization(for test purpose)
% if exist('labels','var')
%     visulizeGroups(x,y,labels+1);
%     figure
%     selected=ones(length(x),1);
%     selected(gated_idx)=2;
%     visulizeGroups(x,y,selected);
%     figure
%     visulizeGroups(x,y,labels);
%     plot(boundaryx,boundaryy);
% end