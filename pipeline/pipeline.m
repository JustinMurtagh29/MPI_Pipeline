function pipeline(p,pT,todo);
% Central function to execute certain parts of the pipeline;
% Create parameters by using [p pT] = setParameterSettings();, afterwards pass p and pT here 
% todo can be any of strings parsed here
% All functions called can be found in pipeline subdirectory of manuelCode repositorium as well

% Classification of dense or training regions
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
	% Should be placed in a single subdirectory as well (as it is not really part of the active classifier)
	case 'graphFeatures'
		miniFeature(p);
	case 'graphFeaturesLR'
		miniFeature(pT);
	% These last ones reside in the active repo and implement the GP classifier in an active fashion 
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

