function pipeline(p,pT,todo);
% Central function to execute certain parts of the pipeline;
% Create parameters by using [p pT] = setParameterSettings();, afterwards pass p and pT here 
% todo can be any of strings parsed here

% Classification of dense or training regions
switch todo
	case 'classification'
		bigFwdPass(p);
	case 'classificationLR'
		minicubeFwdPass(pT);
	case 'segmentation'
		miniSegmentation(p);
	case 'segmentationLR'
		miniSegmentation(pT);
	case 'graphConstruction'
		graphConstruction(p);
	case 'graphConstructionLR'
		graphConstruction(pT);
	case 'correspondence'
		correspondenceFinder(p);
	case 'graphFeatures'
		% Alex please enter here
	case 'prepareSupervoxelGP'
		prepareTrainingData(pT);
		prepareHyperparameter(pT);
	case 'applySupervoxelGP'
		fromGraphToDB(p);
	case 'updateSupervoxelGP'

	otherwise 
		display('Unknown instructions! No actions performed ...');
end

end

