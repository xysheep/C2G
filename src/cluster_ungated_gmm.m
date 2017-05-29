function newlabel = cluster_ungated_gmm(data, label, maxn, minn)
if sum(label==0) == 0
    newlabel = label;
    return
end
n = ceil(sum(unique(label)>0) * length(label) / sum(label>0));
if exist('maxn','var')
    if n > maxn
        n = maxn;
    end
end
if exist('minn','var')
    if n < minn
        n = minn;
    end
end
newlabel = label;

%gmmfit = fitgmdist(data(label==0,:),n);
%newlabel(label==0) = max(label) + cluster(gmmfit,data(label==0,:));
gmmfit = kmeans(data(label==0,:),n);
newlabel(label==0) = max(label) + gmmfit;