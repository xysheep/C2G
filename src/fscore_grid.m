function score = fscore_grid(numc,target,gated_idx)
tp = sum(sum(numc(gated_idx,target+1)));
fp = sum(sum(numc(gated_idx,:))) - tp;
fn = sum(sum(numc(:,target+1))) - tp;
score = 2*tp/(2*tp+fp+fn);
end