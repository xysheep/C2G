function [new_labels,f] = cluster_ungated(data,l,varargin)
%% Function description
% This function is used to cluster ungated cells based their overlap with
% kown target populations. The returned cells labeled 0 in original labels
% are changed into some cluster labels larger than original labels. 
%       new_l = cluster_ungated(data,l,...) "d" is the M-by-N data matrix where M is the
%       number of markers. "l" are M-by-1 matrix represent cell labels before pre-cluster. 
% Optional parameter:
%       'pre_cluster_perc' Minimum percentage of unlabeled cells that
%       should be pre-clustered. Default is 0.95. 
%       'outlier_level'Percentage of low density cells of target
%       populations ignored when decide whether a unlabeled cells overlaped
%       with target populations. Default is 0.05.

%% Address input parameters
if sum(l==0)==0
    new_labels = l;
    return
end
pnames = {'maximum_cluster', 'ignore_small_bin','cluster_amount','outliers','visualize'};
dflts  = {300, 20, 0.95, 0.05, false};
[maximum_cluster,small_bin, pre_cluster_perc,outliers_level,visualize] = internal.stats.parseArgs(pnames,dflts,varargin{:});
%% Compute overlap for target population on each dimension
% One cell is defined to be overlapped with one target population on
% certain dimension if its intensity on that dimension is between maximum
% and minimum intensity of target population on that dimension
unique_l = unique(l);
unique_l(unique_l==0) = [];
n_cluster = length(unique_l);
in_range = zeros(size(data,1),n_cluster*size(data,2));
outliers_level = 1 - outliers_level;
if visualize == 1
    f=figure('Position',[680 678 350 300]);
end
k=1;
for i = 1:length(unique_l)
    for d = 1:size(data,2)
        
        
        %in_range(:,i*(n_cluster-1)+d) = data(:,d)<=max(data(l==i,d)) &...
        %    data(:,d)>=min(data(l==i,d));unique_l(i)
        right = prctile(data(l==unique_l(i),d),...
            100-outliers_level/2);
        left = prctile(data(l==unique_l(i),d),...
            outliers_level/2);
        in_range(data(:,d)>right,(i-1)*size(data,2)+d) = 1;
        in_range(data(:,d)<left,(i-1)*size(data,2)+d) = -1;
        
        if visualize == 1
            subplot(size(data,2),length(unique_l),(d-1)*length(unique_l)+i);
            ksdensity(data(l==unique_l(i),d));
            hold on
            ksdensity(data(l==0,d));
            plot([left left],[0 0.6],'r-')
            plot([right right],[0 0.6],'r-')
            text((left-5)/2, 0.5,'-1');
            text((right+15)/2, 0.5,'+1');
            text((left+right)/2, 0.5,'0');
            %legend('Target','Unknown')
            title(sprintf('Target %d in Marker %d',unique_l(i),d),'FontSize', 12);
            axis([-5 15 0 0.6]);
        end
        k = k + 1;
    end
end

%% Decide number of clusters (remove low frequency pattern)
% Here, we do not really want the exactly number of ungated populations. We
% just want to label out main ungated population, so we can identify manual
% gating sequence easier. For this reason, pick up a hard threshold is OK
% currently. By default, patterns with cummulutive frequency less than 0.1 are ignored
% and considered outliers.

ungated_in_range = in_range(l==0,:); % choose all ungated cells
[~,~,ic] = unique(ungated_in_range,'rows');
unique_pattern = unique(ic);
freq_pattern = histc(ic,unique_pattern);
sorted_freq_pattern = sort(freq_pattern,1,'descend');
total_cells = size(ungated_in_range,1);
s = 0;
for i = 1:length(sorted_freq_pattern)
    s = s + sorted_freq_pattern(i);
    if s >= total_cells * pre_cluster_perc || i > maximum_cluster || sorted_freq_pattern(i)<small_bin%(length(unique(l))-1)*size(data,2)
        break
    end
end
used_patterns = unique_pattern(freq_pattern >= sorted_freq_pattern(i));
new_labels = l;
ungated_labels = zeros(sum(l==0),1);
last_label = max(l) + 1;
for pt = used_patterns'
    ungated_labels(ic==pt) = last_label;
    last_label = last_label +1;
end
new_labels(l==0) = ungated_labels;