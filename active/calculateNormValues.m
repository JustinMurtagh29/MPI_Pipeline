
function calculateNormValues(p,i,j,k)

normFile = [p.saveFolder 'state/normValues.mat'];
load(p.local(i,j,k).weightFile);
X=weights;
percentiles = [.1 .9];
X_lower_cutoff = quantile(X,percentiles(1),1);
X_upper_cutoff = quantile(X,percentiles(2),1);

save(normFile, 'X_lower_cutoff', 'X_upper_cutoff');

end
