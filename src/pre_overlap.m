function [overlap2d, overlap] = pre_overlap(data, means, covs, label)
n_dim = size(data,2);
overlap2d = zeros(length(means),length(means),n_dim,n_dim);
for i = 1:n_dim
    for j = i+1:n_dim
        adj = gaussian_pop_overlap(means,covs,label,[i j],1);
        overlap2d(:,:,i,j) = adj;
        overlap2d(:,:,j,i) = adj;
    end
end
overlap = zeros(length(means));
for i = 1:length(means)
    for j = i+1:length(means)
        overlap(i,j) = bhattacharyya(means{i},means{j},covs{i},covs{j});
    end
end
overlap = overlap + overlap';
        