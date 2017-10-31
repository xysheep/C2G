function noise_l = make_noise(ori_l,percentage)
%% Function make_noise
% This function is used to make noise from origianl labels. For each target
% populations, this function will randomly assign some ungated cells to it
% as noise. 
if ~exist('percentage','var')
    percentage = 1;
end

uniq_labels = unique(ori_l(ori_l~=0));
noise_l = ori_l;
for i = 1:length(uniq_labels)
    idx_ungated = find(noise_l==0);
    num_noise = round(percentage/100 * sum(noise_l == uniq_labels(i)));
    if num_noise < length(idx_ungated)
        selected_idx = randsample(idx_ungated, num_noise);
    else
        selected_idx = idx_ungated;
    end
    noise_l(selected_idx) = uniq_labels(i);
end


