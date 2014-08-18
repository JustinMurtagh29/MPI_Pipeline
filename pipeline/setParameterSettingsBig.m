function [parameter, parameterTrain] = setParameterSettings(old_datestr)
% Function to set parameter for pipeline, give no argument if you want to create new parameters or a former datestring if you want to overwrite 

% GLOBAL SETTINGS FOR graphConstruction.m
% Start time as unique identifier for reference and storage
if nargin == 0
	parameter.start = datestr(clock, 30);
elseif nargin == 1
	parameter.start = old_datestr;
end
parameter.saveFolder = ['/zdata/manuel/results/pipeline/' parameter.start '/'];
% Define region to put through pipeline
parameter.bbox = [641 8320; 769 5888; 129 3200]; % this should be aligned with KNOSSOS cubes and be divisble by tileSize 
%parameter.bbox = [3073 5120; 3073 5120; 2049 3072];
parameter.tileSize =  [512; 512; 256]; % Size of local segmentation and local graph construction
parameter.tileBorder = [-256 256; -256 256; -128 128]; % border of local segmentation included for gloablization and large size due to games
parameter.tiles = (parameter.bbox(:,2) - parameter.bbox(:,1) + 1) ./ parameter.tileSize;
% Which raw dataset
parameter.raw.root = '/zdata/manuel/data/cortex/2012-09-28_ex145_07x2_corrected/color/1/';
parameter.raw.prefix = '2012-09-28_ex145_07x2_corrected_mag1';
% Which classifier to use
parameter.cnn.dateStrings = '20130516T204040';
parameter.cnn.iter = 8; 
parameter.cnn.gpu = 3;
parameter.cnn.first = ['/zdata/manuel/results/parameterSearch/' parameter.cnn.dateStrings '/iter' num2str(parameter.cnn.iter, '%.2i') '/gpu' num2str(parameter.cnn.gpu, '%.2i') '/saveNet0000000001.mat'];
parameter.cnn.GPU = true;
% Function to use for classification
parameter.class.func = @bigFwdPass;
% Location to store CNN classification
parameter.class.root = [parameter.saveFolder 'class/'];
parameter.class.prefix = parameter.raw.prefix;
% Function to use for segmentation
parameter.seg.func = @seg20130604;
% Specify arguments for filterbank applied to raw and aff data each
parameter.filter = {{'sortedeigenvalueshessian' [3 5] []}...
	{'gaussiansmoothedgradmagnitude' [3 5] []}...
	{'intensitygaussiansmoothed' [3 5] []}...
	{'sortedeigenvaluesstructure' [3 5] [5 7]}...
	{'laplaceofgaussian' [3 5] []}...
	{'differenceofgaussians' [3 5] []}};

% Correspondence parameter
parameter.correspondence.overlap = 1; % overlap of local segmentation to compare on each side around a face
parameter.correspondence.saveFolder = [parameter.saveFolder 'correspondences/'];

% Feature parameter
parameter.feature.root = [parameter.saveFolder 'features\'];
% Function to use for FeatureCalculation
parameter.feature.func = @calcFeatures;
% Choice of filters for 'raw' or 'class'
parameter.feature.input = {'raw', 'aff'};

% GLOBAL SETTINGS FOR fromGraphToDB.m
% State variables from the GP
parameter.gp.stateFolder = [parameter.saveFolder 'state/'];
parameter.gp.normValues = [parameter.gp.stateFolder 'normValues.mat'];
parameter.gp.hyperParameter = [parameter.gp.stateFolder 'hyperParameter.mat'];
parameter.gp.initalGroundTruth = [parameter.gp.stateFolder 'initalGroundTruth.mat'];
% Define cutoff(s) for writing to knowledge DB 
parameter.gp.upperCut = .95;
parameter.gp.lowerCut = .15;

% GLOBAL SETTINGS FOR writeKnowledgeDB.m
% Where to put knowledge DB data
parameter.kdb.writeLocation = [parameter.saveFolder 'kdb/'];
parameter.kdb.counter = [parameter.saveFolder 'kdb/counter.mat'];
% Settings for dataset (written to knowledge DB)
parameter.kdb.settings.name = '2012-09-28_ex145_07x2_corrected';
parameter.kdb.settings.scale = [11.24 11.24 28];
parameter.kdb.settings.priority = 0;
% Viewports
% First for TRacer
parameter.kdb.v(1).name = 'tracer';
parameter.kdb.v(1).width = 200;
parameter.kdb.v(1).height = 200;
parameter.kdb.v(1).isRotated = logical(0);
parameter.kdb.v(1).openingAngle = [];
% Second for B4B
parameter.kdb.v(2).name = 'b4b';
parameter.kdb.v(2).width = 312;
parameter.kdb.v(2).height = 214;
parameter.kdb.v(2).isRotated = logical(1);
parameter.kdb.v(2).openingAngle = 20;
%kdb flags
parameter.kdb.saveForProblemInspector = 1;

% LOCAL SETTINGS for each tile
for i=1:parameter.tiles(1)
	for j=1:parameter.tiles(2)
		for k=1:parameter.tiles(3)
			% Save path for data relating to this tile
			parameter.local(i,j,k).saveFolder = [parameter.saveFolder 'x' num2str(i, '%.4i') 'y' num2str(j, '%.4i') 'z' num2str(k, '%.4i') '/'];
			% Bounding box without and with border for this tile
			parameter.local(i,j,k).bboxSmall = [parameter.bbox(:,1) + [i-1; j-1; k-1] .* parameter.tileSize parameter.bbox(:,1) + [i; j; k] .* parameter.tileSize - [1; 1; 1]];
			parameter.local(i,j,k).bboxBig = parameter.local(i,j,k).bboxSmall + parameter.tileBorder;
			% Where to save
			parameter.local(i,j,k).segFile = [parameter.local(i,j,k).saveFolder 'seg.mat'];
			parameter.local(i,j,k).edgeFile = [parameter.local(i,j,k).saveFolder 'edges.mat'];
			parameter.local(i,j,k).borderFile =  [parameter.local(i,j,k).saveFolder 'borders.mat'];
			parameter.local(i,j,k).weightFile = [parameter.local(i,j,k).saveFolder 'weights.mat'];
			parameter.local(i,j,k).probFile = [parameter.local(i,j,k).saveFolder 'prob.mat'];
		end
	end
end


% GLOBAL SETTINGS FOR training data generation
parameterTrain = parameter;
parameterTrain = rmfield(parameterTrain, 'local');
parameterTrain.cnn.GPU = false; % minicubeFwdPass on CPU to allow arbitrary region size

% Regions needed for training data

% Region from Heiko, training region
parameterTrain.local(1).bboxSmall = [4097 4736; 4481 5248; 2250 2450]; % Heiko
parameterTrain.local(1).trainFile = '/zdata/manuel/data/cortex/denseSkel/region1.nml';
parameterTrain.local(1).useSegTraining = true;
parameterTrain.local(1).useGpTraining = true;
% Region from Alex, training region
parameterTrain.local(2).bboxSmall = [1417 1717; 4739 5039; 890 1190]; % Alex
parameterTrain.local(2).trainFile = '/zdata/manuel/data/cortex/denseSkel/region2.nml';
parameterTrain.local(2).useSegTraining = true;
parameterTrain.local(2).useGpTraining = true;
% Region from Max & Anna, test region
parameterTrain.local(3).bboxSmall = [6800 7100; 2140 2440; 1236 1536]; % Max/Anna
parameterTrain.local(3).trainFile = '/zdata/manuel/data/cortex/denseSkel/region3.nml';
parameterTrain.local(3).useSegTraining = false;
parameterTrain.local(3).useGpTraining = false;

% LOCAL SETTINGS FOR training tiles
for i=1:3
	% Save path for data relating to this tile
	parameterTrain.local(i).saveFolder = [parameterTrain.saveFolder 'train' num2str(i, '%.4i') '/'];
        % Bounding box without and with border for this tile
        parameterTrain.local(i).bboxBig = parameterTrain.local(i).bboxSmall + parameterTrain.tileBorder;
        % Where to save
	parameterTrain.local(i).class.root = [parameterTrain.local(i).saveFolder 'class/'];
	parameterTrain.local(i).class.prefix = parameterTrain.class.prefix;
	parameterTrain.local(i).segFile = [parameterTrain.local(i).saveFolder 'seg.mat'];
	parameterTrain.local(i).edgeFile = [parameterTrain.local(i).saveFolder 'edges.mat'];
        parameterTrain.local(i).borderFile =  [parameterTrain.local(i).saveFolder 'borders.mat'];
        parameterTrain.local(i).weightFile = [parameterTrain.local(i).saveFolder 'weights.mat'];
	parameterTrain.local(i).gtFile = [parameterTrain.local(i).saveFolder 'region' num2str(i) '.mat'];
end

% Create state folder for this parameter settings GP
if ~exist(parameter.gp.stateFolder, 'dir')
	mkdir(parameter.gp.stateFolder);
end
% Create folder for writing Knowledge DB for games
if ~exist(parameter.gp.stateFolder, 'dir')
	mkdir(parameter.gp.stateFolder);
end

% Save everything
pT = parameterTrain;
p = parameter;
save([parameter.saveFolder 'allParameter.mat'], 'p', 'pT');

end

