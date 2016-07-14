function makeInterfacePredictions(p)

%Calculate feature quantiles
m=load('/gaba/u/bstaffle/data/Classifier/20150708T153215/20150708T153215Pred.mat'); % Benedikt's Classifier
classifier = m.classificationStruct;
featureMap = m.options.featureMap;
areaT = featureMap.areaThreshold;

for i=1:numel(p.local)

	%load interface wieghts
	m=load([p.local(i).saveFolder 'interfaceWeights.mat']);
	X=m.X;

	% Normalize the features to map the quantiles onto 07x2 feature quantiles
	X= normalizeDataForIC(X,p);

	% Find intIdx i.e. the borders that were considered for interface calculation
	m=load(p.local(i).borderFile);
	borders = m.borders;
	area = [borders(:).Area];
	intIdx = area > areaT;	

	%make predictions
	scores = NaN(length(intIdx),2);
	[~,intScores] = predict(classifier,X);
	scores(intIdx,:) = reshape(intScores,length(intScores)/2,2);

	%Save scores to sunpases.mat file
	m = matfile([p.local(i).saveFolder 'synapses.mat'], 'Writable', true);
	m.scores = scores;

	clear m area borders intIdx scores X 
end

end
