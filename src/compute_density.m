function density = compute_density(data,label)
n_dim = size(data,2);
density = zeros(size(data,1),n_dim,n_dim);
all_labels = unique(label);
all_labels(all_labels==0) = [];

for i = 1:n_dim-1
    for j = i+1:n_dim
%         fprintf('Compute Local-Density for marker %d and marker %d\n',i,j);
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