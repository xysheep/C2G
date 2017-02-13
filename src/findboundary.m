function [boundaryx,boundaryy] = findboundary(ori_xg,ori_yg,x,y)
%% Function findboundary
% This function is used to estimate boundary of a set of data in 2D
% Compared to directly compute convex hull, this function consider the
% case where the input data are colinear. In addition, this function try to
% reduce the number of node in the boundaryx using DouglasPeucker algorithm.


%%  Scripts to test this function
%%
%   Use the following scripts to test this function:
%   
%   fdname = '../FlowData/Primary Murine T cell Data';
%   [data,labels]=load_mul_fcs(fdname,'ctr.fcs');
%   i=8;j=27;
%   figure;colormap(gcf,jet);
%   scatplot(data(:,i),data(:,j),'voronoi',[],100,5,2,4);
%   hold on
%   [bx,by] = findboundary(data(labels==2,i),data(labels==2,j));
%   plot(bx,by,'LineWidth',2);


uniq_d = unique([ori_xg ori_yg],'rows');
xg = uniq_d(:,1);
yg = uniq_d(:,2);
num_uniq = length(xg);

%% Main part
if num_uniq == 1
    [~,d] = knnsearch([x y],[xg yg]);
    boundaryx = [xg+d xg xg-d xg];
    boundaryy = [yg yg-d yg yg+d];
elseif num_uniq >= 3  && rank([xg yg])>1
    k = convhull(xg,yg);
    d = DouglasPeucker([xg(k),yg(k)], 'AUTO', 0);
    [boundaryx,boundaryy] = poly2cw(d(:,1),d(:,2));
elseif num_uniq >= 3 || num_uniq ==2
    [~,beg_idx] = min(xg);
    [~,end_idx] = max(xg);
    b = [xg(beg_idx) yg(beg_idx)];
    e = [xg(end_idx) yg(end_idx)];
    mid = (b+e)/2+(b-e)/10;
    boundaryx = [b(1) e(1) mid(1)];
    boundaryy = [b(2) e(2) mid(2)];
end   
[boundaryx,boundaryy] = poly2cw(boundaryx,boundaryy);

%% Further improvement
% When the input data is not in a perfect convex shape, this function won't
% give a good boundary (may include a large empty region). To address this
% issue, one improvement I can do is when number of input points is small,
% continue to use convex hull. When the numnber of points is large, use
% alternative approach like nearest neighbour or grids.