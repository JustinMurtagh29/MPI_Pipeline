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
    % Segmentation parameter search related
    case 'prepareSkeletons'
        prepareSkeletons(pT);
        skeletonStatistics('/zdata/manuel/data/cortex/denseSkel/');
    case 'segmentationPS'
        parameterSearchSeg(pT);
 	% See graphConstruction subdirectory
	case 'graphConstruction'
		graphConstruction(p);
	case 'graphConstructionLR'
		graphConstruction(pT);
	% See filterbank subdirectory
	case 'filterbank'
		miniFeature(p);
	case 'filterbankLR'
		miniFeature(pT);
    % Globalization of segmentation related    
	case 'correspondence'
		correspondenceFinder(p);
    case 'globalID'
        calculateGlobalID(p);
	% FROM HERE: GP/Supervoxel graph related
	case 'prepareTrainingData'
        prepareTrainingData(pT,'edges');
	% These last ones reside in the active repo and implement the GP classifier in an active fashion
	case 'prepareGP'
		optimizeHyperparameterGP(pT,'edges');
	case 'applyGP'
		makePredictions(p,'edges');
	case 'constructSupervoxelGraph'
		constructSupervoxelGraph(p);
    % Glia prediction using GP functions    
    case 'addNeighborFeatures'
        addNeighborFeatures(p)
    case 'addNeighborFeaturesLR'
        addNeighborFeatures(pT)
    case 'prepareTrainingData Glia'
        prepareTrainingData(pT,'glia');
	case 'prepareGP Glia'
		optimizeHyperparameterGP(pT,'glia');
	case 'applyGP Glia'
		makePredictions(p,'glia');
	otherwise 
		display('Unknown instructions! No actions performed ...');
end

end

