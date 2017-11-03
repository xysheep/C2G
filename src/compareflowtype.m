function bestfscore = compareflowtype(data, label)
threshold = [0.8867975 3.1624543 2.0743298 2.8841788 1.4826645 1.7816912 0.8783884 0.8651572 2.0914003 0.9573006 1.8953659 1.5105588];
lft = data > threshold;
uni_l = unique(label);
uni_l(uni_l == 0) = [];
gidxs = cell(0);
for j = 1:length(uni_l)
    gidxs{j} = (label == uni_l(j));
end
bestfscore = zeros(size(uni_l));
bestgate = zeros(size(uni_l, 1), 2);
t = 0;
for i = 1:(2^size(data,2)-1)
    idx = (dec2bin(i) - '0') == 1;
    idx = [false(1,size(data,2)-length(idx)) idx];
    for k = 0:(2^sum(idx)-1)
        t = t + 1;
        gate = (dec2bin(k) - '0');
        gate = [zeros(1,sum(idx)-length(gate)) gate] == 1;
        gateidx = ismember(lft(:,idx), gate, 'rows');
        if mod(t,1000) ==0
            fprintf('i=%4d\tk=%4d\t%d\n',i,k,sum(gateidx));
        end
        if sum(gateidx) < 10
            continue
        end
        for j = 1:length(uni_l)
            groupidx = gidxs{j};
            tp = sum(groupidx & gateidx);
            fp = sum(~groupidx & gateidx);
            fn = sum(groupidx & ~gateidx);
            fscore = 2 * tp /(2 * tp + fp + fn);
            if fscore > bestfscore(j)
                bestfscore(j) = fscore;
                bestgate(j,:) = [i,k];
            end
        end
    end
end