% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
outDir = '/tmpscratch/amotta/l4/2019-10-22-synem-training-and-test-data-as-hdf5';

% NOTE(amotta): https://gitlab.mpcdf.mpg.de/connectomics/Benedikt/blob/ce5d1b634bb347d953fdc3df7492dfbb2eb63012/+SynEM/+Training/evalPipelineClassifier_v2.m#L47
synemRootDir = '/gaba/u/bstaffle/data/SynEM/data';

% NOTE(amotta): /mnt/mpibr/data/Data is //storage.hest.corp.brain.mpg.de/data
thesisRootDir = '/mnt/mpibr/data/Data/stafflerb/fromPersonal/data/SynEMv2/TestSet';
thesisDirs = {'Cube1082', 'Cube1579'};

% NOTE(amotta): The EM data we've provided with the training data has a
% padding of (100, 100, 40) in both the positive and negative direction
% around the box of interest. See +SynEM/+Script/exportSynapsesAsHdf5.m
paddingBox = [100, 100, 40];

thesisOutPrefix = '2012-09-28_ex145_07x2_ROI2017_TestSet';

l4Roi2017RootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';

info = Util.runInfo();
Util.showRunInfo(info);

%% Test set generated by Benedikt Staffler for the SynEM paper
clear cur*;

% See https://gitlab.mpcdf.mpg.de/connectomics/Benedikt/blob/ce5d1b634bb347d953fdc3df7492dfbb2eb63012/+SynEM/+Training/evalPipelineClassifier_v2.m#L47
curLabels = SynEM.TestSet.loadInterfaceGT(synemRootDir, true, false, 'v3');

% See https://gitlab.mpcdf.mpg.de/connectomics/Benedikt/blob/ce5d1b634bb347d953fdc3df7492dfbb2eb63012/+SynEM/+Training/evalPipelineClassifier_v2.m#L55
curLabels.interfaceLabels = uint32(curLabels.group(:));
curLabels = rmfield(curLabels, 'group');

curSynapseIds = uint32(curLabels.interfaceLabels(:));
curIsInhibitory = uint8(curLabels.isInh);

% See https://gitlab.mpcdf.mpg.de/connectomics/Benedikt/blob/ce5d1b634bb347d953fdc3df7492dfbb2eb63012/+SynEM/+TestSet/loadInterfaceGT.m#L59
curParam = SynEM.TestSet.loadSegParams(synemRootDir);
curCubeParam = curParam.local(67);

curBox = curCubeParam.bboxSmall;
curBoxPadded = curBox + [-1, +1] .* paddingBox(:);

curRaw = loadRawData(curParam.raw, curBoxPadded);
curSeg = loadSegDataGlobal(curParam.seg, curBoxPadded);

curEdges = Util.load(curCubeParam.edgeFile, 'edges');
curEdges = curEdges(curLabels.localEdgeIdx, :);
assert(isequal(size(curEdges, 1), numel(curLabels.interfaceLabels)));
assert(issortedrows(curEdges));

curBorders = reshape({curLabels.borders.PixelIdxList}, [], 1);
curBorders = padBorders(curBox, curBoxPadded, curBorders);
curBorders = cellfun(@uint32, curBorders, 'UniformOutput', false);

% Build output
curOut = struct;
curOut.em = curRaw;
curOut.segmentation = curSeg;
curOut.interfaces = curBorders;
curOut.edges = curEdges;
curOut.synapseIds = curSynapseIds;
curOut.isInhibitory = curIsInhibitory;
    
% Uncomment for debugging
% showMovie(curOut);

% Save output
curOutFile = '2012-09-28_ex145_07x2_segNew_TestSet_Cube67.hdf5';
curOutFile = fullfile(outDir, curOutFile);

structToHdf5(curOutFile, '/', curOut);
infoToHdf5(curOutFile, info);
Util.protect(curOutFile);

%% Test set generated by Benedikt Staffler for his PhD thesis
clear cur*;

curParam = fullfile(l4Roi2017RootDir, 'allParameter.mat');
curParam = Util.load(curParam, 'p');

for curThesisDirIdx = 1:numel(thesisDirs)
    curThesisDir = thesisDirs{curThesisDirIdx};
    curThesisDir = fullfile(thesisRootDir, curThesisDir);
    
    curLabels = load(fullfile(curThesisDir, 'labels.mat'));
    curSynapseIds = uint32(curLabels.interfaceLabels(:));
    
    curBox = curLabels.bbox;
    curBoxPadded = curBox + [-1, +1] .* paddingBox(:);
    
    curRaw = loadRawData(curParam.raw, curBoxPadded);
    curSeg = loadSegDataGlobal(curParam.seg, curBoxPadded);
    
    curCubeIds = curBox(:, 1) - curParam.bbox(:, 1) + 1;
    curCubeIds = ceil(curCubeIds ./ curParam.tileSize);
    
    curCubeId = num2cell(curCubeIds);
    curCubeId = sub2ind(size(curParam.local), curCubeId{:});
    curCubeParam = curParam.local(curCubeId);
    
    curEdges = Util.load(curCubeParam.edgeFile, 'edges');
    assert(isequal(size(curEdges, 1), numel(curLabels.interfaceLabels)));
    assert(issortedrows(curEdges));
    
    curBorders = Util.load(curCubeParam.borderFile, 'borders');
    assert(isequal(numel(curBorders), numel(curLabels.interfaceLabels)));
    
    % Extract voxel indices
    curBorders = reshape({curBorders.PixelIdxList}, [], 1);
    curBorders = padBorders(curBox, curBoxPadded, curBorders);
    curBorders = cellfun(@uint32, curBorders, 'UniformOutput', false);
    
    % Build output
    curOut = struct;
    curOut.em = curRaw;
    curOut.segmentation = curSeg;
    curOut.interfaces = curBorders;
    curOut.edges = curEdges;
    curOut.synapseIds = curSynapseIds;
    
    % Uncomment for debugging
    % showMovie(curOut);
    
    % Save output
   [~, curOutFile] = fileparts(curThesisDir);
    curOutFile = sprintf('%s_%s.hdf5', thesisOutPrefix, curOutFile);
    curOutFile = fullfile(outDir, curOutFile);
    
    structToHdf5(curOutFile, '/', curOut);
    infoToHdf5(curOutFile, info);
    Util.protect(curOutFile);
end

%% Utilities
function borders = padBorders(boxOld, boxNew, borders)
    sizOld = 1 + transpose(diff(boxOld, 1, 2));
    sizNew = 1 + transpose(diff(boxNew, 1, 2));
    padding = boxOld(:, 1) - boxNew(:, 1);
    assert(all(padding >= 0));
    
    for curIdx = 1:numel(borders)
        curIds = {borders{curIdx}(:), [], []};
       [curIds{:}] = ind2sub(sizOld, curIds{1});
       
        for i = 1:3; curIds{i} = curIds{i} + padding(i); end
        
        curIds = sub2ind(sizNew, curIds{:});
        curIds = reshape(curIds, [], 1);
        borders{curIdx} = curIds;
    end
end

function showMovie(out)
    mov = reshape(out.em, [1, size(out.em)]);
    mov = reshape(repmat(mov, 3, 1, 1, 1), 3, []);
    
    curSynVoxels = out.interfaces(out.synapseIds > 0);
    curSynVoxels = cell2mat(curSynVoxels);
    
    % Set synaptic interfaces to red
    mov(1, curSynVoxels) = intmax(class(mov));
    mov(2, curSynVoxels) = 0;
    mov(3, curSynVoxels) = 0;
    
    mov = reshape(mov, [3, size(out.em)]);
    mov = permute(mov, [2, 3, 1, 4]);
    implay(mov);
end