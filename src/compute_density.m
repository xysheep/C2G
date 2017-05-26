function [means, covs, density] = compute_density(data,label)
n_dim = size(data,2);
all_labels = unique(label);
all_labels(all_labels==0) = [];
means = cell(0);
covs = cell(0);
density = zeros(size(data,1), n_dim, n_dim);
% 
% for i = 1:n_dim-1
%     for j = i+1:n_dim
% %         fprintf('Compute Local-Density for marker %d and marker %d\n',i,j);
% %         tic
%         for t = 1:length(all_labels)
%             idx = (label == all_labels(t));
%             dd = zeros(sum(idx),1);
%             rng(9464)
%             x = data(idx,i) + rand(sum(idx),1)*1e-10;
%             y = data(idx,j) + rand(sum(idx),1)*1e-10;
%             if length(x) > 10
%                 [v,c] = voronoin([x,y]);
%                 for k = 1:length(c)
%                     if all(c{k}>1)
%                         a = polyarea(v(c{k},1),v(c{k},2));
%                         dd(k) = 1/a;
%                     end
%                 end
%                 dp = dd/sum(dd);
%             else
%                 dp = ones(size(x));
%             end
%             density(idx,i,j) = dp;
%             density(idx,j,i) = dp;
%         end
% %         toc
%     end
% end
for l = all_labels'
    means{l} = mean(data(label==l,:));
    covs{l}  = cov(data(label==l,:));
    for i = 1:n_dim
        for j = i+1:n_dim
            density(label==l,i,j) = mvnpdf(data(label==l,:),means{l},covs{l});
            density(label==l,j,i) = density(label==l,i,j);
        end
    end
end

