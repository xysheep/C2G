function [down_data,down_label] = down_sample_perpop(data,label,target_num,shrink_factor)
uniq_label = unique(label);
num_label = histc(label,uniq_label);
p_num_label = (num_label/sum(num_label)).^ shrink_factor;
norm_p_num_label = p_num_label/sum(p_num_label) ;
target_num_label = ceil(target_num * norm_p_num_label);


tmp_data = cell(length(uniq_label),1);
tmp_label = cell(length(uniq_label),1);
for i = 1:length(uniq_label)
    idx = find(label==uniq_label(i));
    if target_num_label(i) < length(idx)
        sampled_idx = randsample(idx,target_num_label(i));
    else
        sampled_idx = idx;
    end
    tmp_data{i} = data(sampled_idx,:);
    tmp_label{i} = label(sampled_idx);
end
down_data = cell2mat(tmp_data);
down_label = cell2mat(tmp_label);



% Test the performance
% figure
% uniq_label = unique(label);
% num_label = histc(label,uniq_label);
% plot(num_label/sum(num_label),'LineWidth',4);
% hold on
% xlabel('Population');
% ylabel('Percentage');
% for shrink_factor = [0 0.1 0.2 0.4 0.5 0.6 0.8 0.9]
%     p_num_label = (num_label/sum(num_label)).^ shrink_factor;
%     norm_p_num_label = p_num_label/sum(p_num_label) ;
%     plot(norm_p_num_label)
% end
% legend(cellstr(num2str([1 0 0.1 0.2 0.4 0.5 0.6 0.8 0.9]')));
% title('Down-sampling distribution under difference shrinkage factor');