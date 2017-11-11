function [boundaryx,boundaryy] = findboundary_grid(ori_xg,ori_yg)
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
% x = [ori_xg;ori_xg+1;ori_xg;ori_xg+1];
% y = [ori_yg;ori_yg+1;ori_yg+1;ori_yg];
% 
% uniq_d = unique([x y],'rows');
uniq_d = unique([ori_xg ori_yg],'rows');
xg = uniq_d(:,1);
yg = uniq_d(:,2);
num_uniq = length(xg);
%% Main part
% Here is a potential mistake. For grid, one grid one gate is allowed.
if num_uniq == 1
    d = 0.25;
    boundaryx = [xg+d xg xg-d xg]';
    boundaryy = [yg yg-d yg yg+d]';
elseif ~colinear(xg,yg)
    k = convhull(xg,yg);
    %d = DouglasPeucker([xg(k),yg(k)], 'AUTO', 0);
    [boundaryx,boundaryy] = poly2cw(xg(k),yg(k));
%     boundaryx = xg(k);
%     boundaryy = yg(k);
else
    [~,beg_idx] = min(xg);
    [~,end_idx] = max(xg);
    b = [xg(beg_idx) yg(beg_idx)];
    e = [xg(end_idx) yg(end_idx)];
    boundaryx = [b(1)+0.5 b(1)-0.5 e(1)-0.5 e(1)+0.5]';
    boundaryy = [b(2) b(2) e(2) e(2)]';
end   
[boundaryx,boundaryy] = poly2cw(boundaryx,boundaryy);

%% Further improvement
% When the input data is not in a perfect convex shape, this function won't
% give a good boundary (may include a large empty region). To address this
% issue, one improvement I can do is when number of input points is small,
% continue to use convex hull. When the numnber of points is large, use
% alternative approach like nearest neighbour or grids.