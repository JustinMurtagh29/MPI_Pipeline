function pipeline(p,pT,todo);
% Central function to execute certain parts of the pipeline;
% Create parameters by using [p pT] = setParameterSettings();, afterwards pass p and pT here 
% todo can be any of strings parsed here
% All functions called can be found in pipeline subdirectory of manuelCode repositorium as well

switch todo
	% Classification (calls fwdPass3DonKnossosFolder from CNN subrepo)
	case 'classification'
		bigFwdPass(p);
	case 'classificationLR'
		minicubeFwdPass(pT);
	% Segmentation (calls function defined in parameter settings file based on segmentation repo)
	case 'segmentation'
		miniSegmentation(p);
	case 'segmentationLR'
		miniSegmentation(pT);
	% See graphConstruction subdirectory
	case 'graphConstruction'
		graphConstruction(p);
	case 'graphConstructionLR'
		graphConstruction(pT);
	% See correspondence subdirectory
	case 'correspondence'
		correspondenceFinder(p);
	% See filterbank subdirectory
	case 'filterbank'
		miniFeature(p);
	case 'filterbankLR'
		miniFeature(pT);
	% FROM HERE: GP/Supervoel graph related
	case 'prepareTrainingData'
		prepareTrainingData(pT);
	% These last ones reside in the active repo and implement the GP classifier in an active fashion
	case 'prepareGP'
		optimizeHyperparameterGP(pT);
	case 'applyGP'
		makePredictions(p);
	case 'constructSupervoxelGraph'
		constructSupervoxelGraph(p);
	otherwise 
		display('Unknown instructions! No actions performed ...');
end

end

