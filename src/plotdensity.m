function plotdensity(x,y,f)
[~,~,rank] = unique(f);
rank = rank/max(rank);
cm = redbluecmap(200); % returns the current color map

colorID = max(1, sum(rank > 0:1/length(cm(:,1)):1,2)); 

myColor = cm(colorID, :); % returns your color
scatter(x,y,10,myColor);