function [best_idx,best_perc,best_boundary] = gated_cells(x,y,gatelabels,target_labels,density)
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
best_idx = 1:length(x);
best_perc = 0;
best_boundary = 0;
for ignore_perc = [0 0.05 0.1 0.2 0.3]
    idx_cell = cell(length(gatelabels),1);
    boundariesx = [];
    boundariesy = [];
    for i = 1:length(gatelabels)
        idx  = gatelabels{i};
        xg = x(idx);
        yg = y(idx);
        [boundaryx,boundaryy] = findboundary_outliers(xg,yg,density(idx),ignore_perc);
        if length(boundaryx) > 2 % At least 3 points to be a polygon
            idx_cell{i} = find(inpolygon(x,y,boundaryx,boundaryy));
            if isempty(boundariesy)
                boundariesx = boundaryx;
                boundariesy = boundaryy;
            else
                [boundariesx,boundariesy] = polybool('union',boundariesx,boundariesy,boundaryx,boundaryy);
            end
        end
    end
    gated_idx = unique(cell2mat(idx_cell));
    
    f_score = fscore(ismember(1:length(x),gated_idx),target_labels');
    if f_score > max_fscore
        max_fscore = f_score;
        best_idx = gated_idx;
        best_perc = ignore_perc;
        best_boundary = [boundariesx boundariesy];
    end
end

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