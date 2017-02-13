function score = fscore(obs,target)
tp = sum(obs & target);
fp = sum(obs & ~target);
fn = sum(~obs & target);
score = 2*tp/(2*tp+fp+fn);
end