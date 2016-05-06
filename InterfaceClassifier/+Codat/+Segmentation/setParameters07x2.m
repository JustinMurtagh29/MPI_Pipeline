function [parameter, parameterTrain] = setParameters07x2(saveFolder, cnnFile, segName)
%SETPARAMETERS Create the folder structure and parameter structs for a
% segmentation pipeline.
% INPUT saveFolder: Full path to folder of the new segmentation. The
%       	segmentation will be saved in a further subfolder (see also
%           segName input). If this folder does not exist it will be
%           created.
%       cnnFile: (Optional) Filename of cnn matfile used for training. The
%       	full path of the file will be in the segmentation main folder
%       	and the matfile has to be copied to the corresponding location.
%       	In case a different location is required this can be changed
%       	manually in the output parameters. The file must contain the
%           saved as variable 'cnet'.
%       	(Default: [parameter.saveFolder 'cnn.mat'])
%           Furthermore there is the cnn.GPU flag used to switch between
%           GPU and CPU computation and the cnn.range field defining the
%           boundary map values (see Segmentation.morphRecon).
%       segName: (Optional) Name for of the segmentation. This can be used
%           to recreate or modify the parameter struct for a previously
%           created segmentation.
%           (Default: datestr(clock,30))
%
% NOTE Folders in the output structs always contain a file separator in the
%      end.
%
% Author: Manuel Berning <manuel.berning@brain.mpg.de>
% Modified by Benedikt Staffler <benedikt.staffler@brain.mpg.de>

if exist('segName','var') && ~isempty(segName)
    parameter.segName = segName;
else %set datestring as unique identifier
    parameter.segName = datestr(clock,30);
end
fprintf('[%s] Creating segmentation %s.\n', datestr(now), parameter.segName);

%set and create segmentation folder
parameter.saveFolder = [Util.addFilesep(saveFolder) parameter.segName filesep];

%set 07x2 related parameters
parameter.raw.root = '/gaba/u/mberning/data/cortex/2012-09-28_ex145_07x2_corrected/color/1/';
parameter.raw.prefix = '2012-09-28_ex145_07x2_corrected_mag1';
parameter.raw.voxelSize = [11.24, 11.24, 28];
parameter.bbox = [641 8320; 769 5376; 129 3200]; % this should be aligned with KNOSSOS cubes and be divisble by tileSize
parameter.tileSize =  [512; 512; 256]; % Size of local segmentation and local graph construction
parameter.tileBorder = [-256 256; -256 256; -128 128]; % border of local segmentation included for gloablization and large size due to games
parameter.tiles = (parameter.bbox(:,2) - parameter.bbox(:,1) + 1) ./ parameter.tileSize;

% membrane classifier
if exist('cnnFile','var') && ~isempty(cnnFile)
    parameter.cnn.saveFile = [parameter.saveFolder cnnFile];
else
    parameter.cnn.saveFile = [parameter.saveFolder 'cnn.mat'];
    warning('CNN save path not specified and will be set to default: %s.',parameter.cnn.saveFile);
end
parameter.cnn.GPU = true;
parameter.cnn.range = [-1, 1];

% Location to store CNN classification
parameter.class.root = [parameter.saveFolder 'class/'];
parameter.class.prefix = parameter.raw.prefix;

% LOCAL SETTINGS for each tile
for i=1:parameter.tiles(1)
    for j=1:parameter.tiles(2)
        for k=1:parameter.tiles(3)
            % Save path for data relating to this tile
            parameter.local(i,j,k).saveFolder = [parameter.saveFolder 'local/' 'x' num2str(i, '%.4i') 'y' num2str(j, '%.4i') 'z' num2str(k, '%.4i') '/'];
            % Bounding box without and with border for this tile
            parameter.local(i,j,k).bboxSmall = [parameter.bbox(:,1) + [i-1; j-1; k-1] .* parameter.tileSize parameter.bbox(:,1) + [i; j; k] .* parameter.tileSize - [1; 1; 1]];
            parameter.local(i,j,k).bboxBig = parameter.local(i,j,k).bboxSmall + parameter.tileBorder;
            % Where to save
            parameter.local(i,j,k).segFile = [parameter.local(i,j,k).saveFolder 'seg.mat'];
            parameter.local(i,j,k).edgeFile = [parameter.local(i,j,k).saveFolder 'edges.mat'];
            parameter.local(i,j,k).borderFile =  [parameter.local(i,j,k).saveFolder 'borders.mat'];
            parameter.local(i,j,k).weightFile = [parameter.local(i,j,k).saveFolder 'weights.mat'];
            parameter.local(i,j,k).probFile = [parameter.local(i,j,k).saveFolder 'prob.mat'];
            % Same files for glia prediction
            parameter.local(i,j,k).segmentFile = [parameter.local(i,j,k).saveFolder 'segments.mat'];
            parameter.local(i,j,k).segmentWeightFile = [parameter.local(i,j,k).saveFolder 'segmentWeights.mat'];
            parameter.local(i,j,k).gliaProbFile = [parameter.local(i,j,k).saveFolder 'gliaProb.mat'];
            %synapse detection
            parameter.local(i,j,k).synapseFile = [parameter.local(i,j,k).saveFolder 'synapses.mat'];
        end
    end
end

% GLOBAL SETTINGS FOR gp training data generation
parameterTrain = parameter;
% Remove all fields that do not make sense in training data setting
parameterTrain = rmfield(parameterTrain, {'local' 'bbox' 'tileSize' 'tiles'});

% Densly skeletonized regions in dataset
% Region from Heiko (SegEM segmentation training)
parameterTrain.local(1).bboxSmall = [4097 4736; 4481 5248; 2250 2450];
% parameterTrain.local(1).trainFileRaw = '/gaba/u/mberning/data/cortex/denseSkel/region1.nml';
% parameterTrain.local(1).trainFileGlia = '/gaba/u/mberning/data/cortex/denseSkel/region1glia.nml';
parameterTrain.local(1).trainFileLocal = '/gaba/u/bstaffle/code/workspace/+Codat/+Segmentation/denseSkel/cortex_training_local_subsampled.nml';
% parameterTrain.local(1).trainFileLocalWithoutGlia = '/gaba/u/mberning/data/cortex/denseSkel/region1localWithoutGlia.nml';

% Region from Alex (SegEM segmentation test)
parameterTrain.local(2).bboxSmall = [1417 1717; 4739 5039; 890 1190];
% parameterTrain.local(2).trainFileRaw = '/gaba/u/mberning/data/cortex/denseSkel/region2.nml';
% parameterTrain.local(2).trainFileGlia = '/gaba/u/mberning/data/cortex/denseSkel/region2glia.nml';
parameterTrain.local(2).trainFileLocal = '/gaba/u/bstaffle/code/workspace/+Codat/+Segmentation/denseSkel/cortex_test_local_subsampled.nml';
% parameterTrain.local(2).trainFileLocalWithoutGlia = '/gaba/u/mberning/data/cortex/denseSkel/region2localWithoutGlia.nml';

% Training region 3 is only required for GP performance evaluation
% % Region from Max & Anna
% parameterTrain.local(3).bboxSmall = [6800 7100; 2140 2440; 1236 1536];
% parameterTrain.local(3).trainFileRaw = '/gaba/u/mberning/data/cortex/denseSkel/region3.nml';
% parameterTrain.local(3).trainFileGlia = '/gaba/u/mberning/data/cortex/denseSkel/region3glia.nml';
% parameterTrain.local(3).trainFileLocal = '/gaba/u/mberning/data/cortex/denseSkel/region3local.nml';
% parameterTrain.local(3).trainFileLocalWithoutGlia = '/gaba/u/mberning/data/cortex/denseSkel/region3localWithoutGlia.nml';

% LOCAL SETTINGS FOR training tiles
for i=1:numel(parameterTrain.local)
    % Save path for data relating to this tile
    parameterTrain.local(i).saveFolder = [parameterTrain.saveFolder 'train' num2str(i, '%.4i') filesep];
    % Bounding box without and with border for this tile
    parameterTrain.local(i).bboxBig = parameterTrain.local(i).bboxSmall + parameterTrain.tileBorder;
    % Where to save
    parameterTrain.local(i).class.root = [parameterTrain.local(i).saveFolder 'class' filesep];
    parameterTrain.local(i).class.prefix = parameterTrain.class.prefix;
    parameterTrain.local(i).seg.parameterSearchFolder = [parameterTrain.local(i).saveFolder 'parameterSearch' filesep];
    parameterTrain.local(i).segFile = [parameterTrain.local(i).saveFolder 'seg.mat'];
    parameterTrain.local(i).edgeFile = [parameterTrain.local(i).saveFolder 'edges.mat'];
    parameterTrain.local(i).borderFile =  [parameterTrain.local(i).saveFolder 'borders.mat'];
    parameterTrain.local(i).weightFile = [parameterTrain.local(i).saveFolder 'weights.mat'];
    parameterTrain.local(i).gtFile = [parameterTrain.local(i).saveFolder 'region' num2str(i) '.mat'];
    
    % Benjamin's glia prediction
    parameterTrain.local(i).segmentFile = [parameterTrain.local(i).saveFolder 'segments.mat'];
    parameterTrain.local(i).segmentWeightFile = [parameterTrain.local(i).saveFolder 'segmentWeights.mat'];
    parameterTrain.local(i).gliaProbFile = [parameterTrain.local(i).saveFolder 'gliaProb.mat'];
end

%create segmentation save folder if it does not exist
if ~exist(parameter.saveFolder, 'dir');
    mkdir(parameter.saveFolder);
    fprintf('[%s] Creating new segmentation folder %s.\n',datestr(now), ...
             parameter.saveFolder);
end

% Save everything
pT = parameterTrain;
p = parameter;
m = matfile([parameter.saveFolder 'allParameter.mat'],'Writable',true);
m.p = p;
m.pT = pT;

end

