function calculateFeatureQuantilesForIC(p)

cubeNo = randperm(numel(p.local),100);
XAll=[];

for i=1:length(cubeNo)
	m=load([p.local(cubeNo(i)).saveFolder 'interfaceWeights.mat']);
	X=m.X;
	XAll = vertcat(XAll,X);
end

%Calculate quantiles 

percentiles = [.1 .9];
X_lower_cutoff = quantile(XAll,percentiles(1),1);
X_upper_cutoff = quantile(XAll,percentiles(2),1);

save([p.saveFolder 'state/featureQuantilesForIC.mat'], 'X_lower_cutoff', 'X_upper_cutoff');

end
