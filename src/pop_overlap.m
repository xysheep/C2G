function [adj,areas,all_labels] = pop_overlap(x,y,label,density,ig_ratio)
all_labels = unique(label);
all_labels(all_labels==0 | histc(label,all_labels)<1) = [];
boundaryx = cell(length(all_labels),1);
boundaryy = cell(length(all_labels),1);
%fprintf('input_c = %d\t',length(x));

if ~exist('mode','var') || strcmp(mode,'polygon')
    idx = {[]};
    for i = 1:length(all_labels)
        idx{i} = (label == all_labels(i));
        xg = x(idx{i});
        yg = y(idx{i});
        [boundaryx{i},boundaryy{i}] = findboundary_outliers(xg,yg,density(idx{i}),ig_ratio);%findboundary(xg,yg,x,y);
    end
    adj = eye(length(all_labels));
    areas = eye(length(all_labels));
    for i = 1:length(all_labels)-1
        for j = i+1:length(all_labels)
            n_int_i = sum(inpolygon(x(idx{i}),y(idx{i}),boundaryx{j},boundaryy{j}));
            n_int_j = sum(inpolygon(x(idx{j}),y(idx{j}),boundaryx{i},boundaryy{i}));
            %adj(j,i) = n_int_i*n_int_j/sum(idx{i})/sum(idx{j});
            adj(j,i) = max(n_int_i/sum(idx{i}),n_int_j/sum(idx{j}));
            areas(j,i) = n_int_i/sum(idx{i});
            areas(i,j) = n_int_j/sum(idx{j});
        end
    end
end
adj = adj + adj' - eye(length(all_labels));
%areas = areas + areas' - eye(length(all_labels));
adj(isnan(adj))=0;
%adj(isnan(areas))=0;

