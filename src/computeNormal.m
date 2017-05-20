function [MUs, SIGMAs] = computeNormal(data, label)
uniquelabel = unique(label);
uniquelabel(uniquelabel == 0) = [];
MUs = cell(max(uniquelabel),1);
SIGMAs = cell(max(uniquelabel),1);
for l = uniquelabel(:)'
    tmpd = data(label == l,:);
    MUs{l} = mean(tmpd);
    SIGMAs{l} = cov(tmpd);
end


