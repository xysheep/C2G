function displayclasify(idx,label,col)
uniq_idx = unique(idx);
size_idx = histc(idx,uniq_idx);
[~,I] = sort(size_idx,2,'descend');
uniq_idx = uniq_idx(I);
% num_balls = length(idx);
% m = zeros(num_balls,num_balls);
% for i = 1:length(uniq_idx)
%     idx_groupi = find(idx==uniq_idx(i));
%     for j = 1:(length(idx_groupi)-1)
%         for k = (j+1):length(idx_groupi)
%             m(idx_groupi(j),idx_groupi(k)) = 1;
%             m(idx_groupi(k),idx_groupi(j)) = 1;
%         end
%     end
% end
% 
% 
% g = graph(m);
% h = plot(g,'layout','subspace');
% set(h,'MarkerSize',5);
% labelnode(h,1:length(label),label');
% set(gca,'xtick',[])
% set(gca,'ytick',[])
% box on;
%axis off;
hold on;

if length(idx)==6
    if max(size_idx)==3
        location = [1 1 1 2 2 3; 1 2 3 1 2 1];
    elseif sum(size_idx==2) == 2
        location = [1 1 2 2 3 3; 1 2 1 2 1 2];
    else
        location = [1 1 3 3 3 3; 1 2 0 1 2 3];
    end
    k = 1;
    for i = 1:length(uniq_idx)
        idx_idx = find(idx == uniq_idx(i));
        for j = 1:length(idx_idx)
            plot(location(1,k),location(2,k),'.',...
                'markersize', 50,'Color', col(label(idx_idx(j))+1,:));
            text(location(1,k),location(2,k),num2str(label(idx_idx(j))),...
                'HorizontalAlignment','center','FontSize',15);
            k = k + 1;
        end
        if length(idx_idx) == 2
            x = location(1,k-1);
            y = location(2,k-1)-0.5;
            drawellipse(x,y,0.4,1);
        elseif length(idx_idx) == 3
            x = location(1,k-1);
            y = location(2,k-1)-1;
            drawellipse(x,y,0.4,1.5);
        else
            x = location(1,k-1);
            y = location(2,k-1);
            drawellipse(x,y,0.4,0.4);
        end
    end
    if max(size_idx)==3
        axis([0 4 0 4]);
    else
        axis([0 4 -0.5 3.5]);
    end
elseif length(idx) == 3
    if max(size_idx) == 2
        location = [1 1 2;1 2 1];
    else
        location = [2 2 2;0.5 1.5 2.5];
    end
    k = 1;
    for i = 1:length(uniq_idx)
        idx_idx = find(idx == uniq_idx(i));
        for j = 1:length(idx_idx)
            plot(location(1,k),location(2,k),'.',...
                'markersize', 50,'Color', col(label(idx_idx(j))+1,:));
            text(location(1,k),location(2,k),num2str(label(idx_idx(j))),...
                'HorizontalAlignment','center','FontSize',15);
            k = k + 1;
        end
        if length(idx_idx) == 2
            x = location(1,k-1);
            y = location(2,k-1)-0.5;
            drawellipse(x,y,0.3,1);
        elseif length(idx_idx) == 3
            x = location(1,k-1);
            y = location(2,k-1)-1;
            drawellipse(x,y,0.3,1.5);
        else
            x = location(1,k-1);
            y = location(2,k-1);
            drawellipse(x,y,0.3,0.3);
        end
    end
    axis([0 3 0 3]);
end
set(gca,'xtick',[])
set(gca,'ytick',[])
box on;