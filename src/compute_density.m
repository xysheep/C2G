function [means, covs, density] = compute_density(data,label)
n_dim = size(data,2);
all_labels = unique(label);
all_labels(all_labels==0) = [];
means = cell(0);
covs = cell(0);
density = zeros(size(data,1), n_dim, n_dim);
% 
fprintf('Compute 2D Local-Density for all marker pairs [');
tmp_n = 0;
for i = 1:n_dim-1
    for j = i+1:n_dim
        tmp_n = tmp_n + 1;
        perc = round(100*2*tmp_n./(n_dim*(n_dim-1)));
        if tmp_n > 1
            fprintf('\b\b\b\b\b')
        end
        fprintf('%3d%%]',perc);
%         tic
        for t = 1:length(all_labels)
            idx = (label == all_labels(t));
            dd = zeros(sum(idx),1);
            rng(9464)
            x = data(idx,i) + rand(sum(idx),1)*1e-10;
            y = data(idx,j) + rand(sum(idx),1)*1e-10;
            if length(x) > 10
                [v,c] = voronoin([x,y]);
                for k = 1:length(c)
                    if all(c{k}>1)
                        a = polyarea(v(c{k},1),v(c{k},2));
                        dd(k) = 1/a;
                    end
                end
                dp = dd/sum(dd);
            else
                dp = ones(size(x));
            end
            density(idx,i,j) = dp;
            density(idx,j,i) = dp;
        end
%         toc
    end
end
fprintf('\n');
% for l = all_labels'
%     means{l} = mean(data(label==l,:));
%     covs{l}  = cov(data(label==l,:));
%     for i = 1:n_dim
%         for j = i+1:n_dim
%             mu = means{l};
%             sigma = covs{l};
%             density(label==l,i,j) = mvnpdf(data(label==l,[i j]),mu([i j]),sigma([i j],[i j]));
%             density(label==l,j,i) = density(label==l,i,j);
%         end
%     end
% end

% for i = 1:n_dim
%     for j = i+1:n_dim
%         bd = cell(length(all_labels),1);
%         for k = 1:length(all_labels)
%             xg = data(label==all_labels(k),i);
%             yg = data(label==all_labels(k),j);
%             [boundaryx,boundaryy] = findboundary(xg,yg);
%             bd{k} = [boundaryx,boundaryy];
%         end
%         for i_g = 1:length(all_labels)-1
%             for j_g = i_g+1:length(all_labels)
%                 if inpolygon(bd{i_g}(:,1),bd{i_g}(:,2),bd{j_g}(:,1),bd{j_g}(:,2))
%                     overlap2d(i_g,j_g,i,j) = true;
%                     overlap2d(j_g,i_g,i,j) = true;
%                     overlap2d(i_g,j_g,j,i) = true;
%                     overlap2d(j_g,i_g,j,i) = true;
%                 end
%             end
%         end
%     end
% end




