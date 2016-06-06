
function calculateNormValues(p)

display('Calculating normValues for GP by sampling weights in random cubes');
normFile = [p.saveFolder 'state/normValues.mat'];

nrCubesToSample = 50;
idxCubesToSample = randi(numel(p.local),nrCubesToSample,1);
X = [];
for i=1:nrCubesToSample
	m=load(p.local(idxCubesToSample(i)).weightFile);
	X =vertcat(X,m.weights);
end

percentiles = [.1 .9];
X_lower_cutoff = quantile(X,percentiles(1),1);
X_upper_cutoff = quantile(X,percentiles(2),1);

save(normFile, 'X_lower_cutoff', 'X_upper_cutoff');

end
