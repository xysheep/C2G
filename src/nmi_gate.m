function [nmi,metric] = nmi_gate(labels,gatelabels)


%% Testing command: 
%% 
% Use the following scripts to test this file 
%   
%   fdname = '../FlowData/Primary Murine T cell Data';
%   [data,labels]=load_mul_fcs(fdname,'ctr.fcs');
%	gatelabels={find(labels==1) find(labels==2) find(labels==3)}
%   x=data(:,i);y=data(:,j);
%	new_entropy_gate(labels,gatelabels)
%	i=8;j=27;
%   gatelabels=new_bestgate(data(:,i),data(:,j),labels);
%   new_entropy_gate(labels,gatelabels)



%% Initiate
% gatelabels = gatelabels(idx);
n = length(labels);
neargroup = false(n,length(gatelabels)+1);

%% Summarize gates into a binary matrix
% Rows correspond to different cells and columns correspond to whether a
% certain cell in certain gates
for i = 1:length(gatelabels)
    neargroup(gatelabels{i},i) = true;
end
neargroup(:,length(gatelabels)+1) = sum(neargroup,2)==0;% If not in any gates, in "ungated" gate

%% Calculate Entropy
norm_neargroup = bsxfun(@rdivide,neargroup,sum(neargroup,2));
% if one cell in mulpitle gates, lower its weight


all_l = unique(labels)';
metric = zeros(length(all_l),length(gatelabels)+1);
% Row of the metric is original cell population(0 is ungated)
% Column of the metric is gates

for il = 1:length(all_l)
        metric(il,:) = sum(norm_neargroup(labels==all_l(il),:),1);
        if (all_l(il)==0)
            metric(il,1:end-1)=metric(il,1:end-1)/2;
            metric(il,end)=metric(il,end)+sum(metric(il,1:end-1),2);
        end
end
%disp(metric)
% entropy_weight = 0;
% n_all = sum(metric(:));
% for i = 1:size(metric,1)
%     weight = sum(metric(i,:));
%     p = metric(i,:)/weight;
%     entropy_weight = entropy_weight - sum(p.*log(p),'omitnan')*weight/n_all;
% end
%p = metric(:)./sum(metric(:));
%entropy = sum(-p.*log(p),'omitnan');


%% Instead of entropy, use NMI
p = metric ./ sum(metric(:));
I = 0;
px = sum(p);
py = sum(p,2);
for i = 1:size(p,1)
    for j = 1:size(p,2)
        if p(i,j) > 0 && px(j) >0 && py(i) >0
            I = I + p(i,j) * log(p(i,j)/px(j)/py(i));
        end
    end
end
if sum(-px.*log(px)) > 0 && sum(-py.*log(py))>0
    nmi = I/sqrt(sum(px.*log(px))*sum(py.*log(py)));
else
    nmi = 0;
end




