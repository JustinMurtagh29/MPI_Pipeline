function X_out = normalizeDataForIC(X)

me=mfilename;
mydir = which(me); mydir = mydir(1:end-2-numel(me));

m=load([mydir '../featureQuantilesForRef.mat']);
refQuantiles = m.quant;

load([p.saveFolder 'state/featureQuantilesForIC.mat']);

X = bsxfun(@minus, X, X_lower_cutoff);
X_norm = bsxfun(@times, X, 1./(X_upper_cutoff-X_lower_cutoff));

X_temp = bsxfun(@times, X_norm,(refQuantiles(2,:) - refQuantiles(1,:)));

X_out = bsxfun(@plus,X_temp,refQuantiles(1,:));

end
