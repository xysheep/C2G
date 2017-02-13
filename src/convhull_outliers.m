function k=convhull_outliers(x,y,outlier)
n = length(x);
idx = 1:n;
while length(idx)>outlier*n
    %fprintf('n=%d outlier=%.2f remains=%d\n',n,outlier,length(idx));
    k = convhull(x(idx),y(idx));
    remain_idx = idx(~ismember(1:length(idx),k));
    if length(unique([x(remain_idx) y(remain_idx)],'rows')) > 3
        idx(k) = [];
    else
        break
    end
end
%fprintf('n=%d outlier=%.2f remains=%d\n',n,outlier,length(idx));
k = convhull(x(idx),y(idx));
k = idx(k);