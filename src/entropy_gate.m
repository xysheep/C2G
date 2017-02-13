function [entropy,num_gates] = entropy_gate(x,y,labels,gatelabels)


%% Testing command: 
%% 
% Use the following scripts to test this file 
%   
%   fdname = '../FlowData/Primary Murine T cell Data';
%   [data,labels]=load_mul_fcs(fdname,'ctr.fcs');
%	i=8;j=27;
%	gatelabels=bestgate(data(:,i),data(:,j),labels,'grids');
%   x=data(:,i);y=data(:,j);
%	entropy_gate(data(:,i),data(:,j),labels,gatelabels)
%



%% Initiate
% idx = (labels~=0);
% x = x(idx);
% y = y(idx);
% labels = labels(idx);
% gatelabels = gatelabels(idx);
n = length(gatelabels);
all_labels = unique(gatelabels)';
neargroup = false(n,length(all_labels));
%% Convex hull
for i = 1:length(all_labels)
    idx = (gatelabels == all_labels(i));
    xg = x(idx);
    yg = y(idx);
    if rank([xg yg])>1
        k = convhull(xg,yg);
        d = DouglasPeucker([xg(k),yg(k)], 'AUTO', 0);
%         boundaryx = d(:,1);
%         boundaryy = d(:,2);
        [boundaryx,boundaryy] = poly2cw(d(:,1),d(:,2));
    else
        [~,beg_idx] = min(xg);
        [~,end_idx] = max(xg);
        b = [xg(beg_idx) yg(beg_idx)];
        e = [xg(end_idx) yg(end_idx)];
        mid = (b+e)/2+(b-e)/10;
        boundaryx = [b(1) e(1) mid(1)];
        boundaryy = [b(2) e(2) mid(2)];
    end  
%     k = convhull(xg,yg);
%     b = DouglasPeucker([xg(k),yg(k)], 'AUTO', 0);
    neargroup(:,i) = inpolygon(x,y,boundaryx,boundaryy);
end
%% Calculate Frequency
% observations = [labels neargroup];
% [~,~,ic] = unique(observations,'rows');
% state_freq = histc(ic,unique(ic));
% p = state_freq/sum(state_freq);



norm_neargroup = bsxfun(@rdivide,neargroup,sum(neargroup,2));
all_l = unique(labels)';
metric = zeros(length(all_l),length(all_labels));
for il = 1:length(all_l)
        metric(il,:) = sum(norm_neargroup(labels==all_l(il),:),1);
end
p = metric(:)./sum(metric(:));
entropy = sum(-p.*log(p),'omitnan');
num_gates = length(unique(gatelabels));
