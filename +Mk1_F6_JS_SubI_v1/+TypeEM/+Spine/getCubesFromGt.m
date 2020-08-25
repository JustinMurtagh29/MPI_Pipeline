% script to extract GT from nmls and their corresponding cubeIds

rootDir = '/tmpscratch/sahilloo/data/Mk1_F6_JS_SubI_v1/pipelineRun_mr2e_wsmrnet';
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;
param.experimentName = 'Mk1_F6_JS_SubI_v1_mrnet_wsmrnet';
maxSegId = Seg.Global.getMaxSegId(param);
import Mk1_F6_JS_SubI_v1.TypeEM.*

addTypeEM = true;

if addTypeEM
    nmlDir = fullfile(param.saveFolder, ...
         'tracings', 'typeEM');
    nmlFiles = fullfile(nmlDir, 'proofread', 'withoutSpines', ...
         {'box-1.nml','box-2.nml', 'box-3.nml', ...
         'box-4.nml','box-5.nml','box-6.nml',...
         'box-7.nml','box-8.nml','box-9.nml',...
         'box-10.nml','box-11.nml','box-12.nml',...
         'box-13.nml'});
    featureSetName = 'segmentAgglomerate'; %'segmentAgglomerate'; % 'segment'
    idxTrain = [1,2,3,4,5,6,7,8,9,10,11,12,13];
    % load features
    Util.log(sprintf('Loading typeEM GT from %s',nmlDir))
    gtType = TypeEM.GroundTruth.loadSet( ...
            param, featureSetName, nmlFiles(idxTrain));

else
    gtType.segId = [];
end

% load SH training data
idxTrain = 1:22;
nmlDir = fullfile(param.saveFolder, ...
     'tracings', 'box-seeded','spine-head-ground-truth', 'v4');
Util.log(sprintf('Loading SH GT from %s', nmlDir))
nmlFiles = fullfile(nmlDir, ...
     {'spine-head-ground-truth-1.nml','spine-head-ground-truth-2.nml', 'spine-head-ground-truth-3.nml', ...
     'spine-head-ground-truth-4.nml','spine-head-ground-truth-5.nml','spine-head-ground-truth-6.nml',...
     'spine-head-ground-truth-7.nml','spine-head-ground-truth-8.nml','spine-head-ground-truth-9.nml',...
     'spine-head-ground-truth-10.nml','spine-head-ground-truth-11.nml',...
     'spine-head-ground-truth-13.nml', ...
     'spine-head-ground-truth-14.nml','spine-head-ground-truth-15.nml','spine-head-ground-truth-16.nml',...
     'spine-head-ground-truth-18.nml','spine-head-ground-truth-19.nml',...
     'spine-head-ground-truth-20.nml', 'spine-head-ground-truth-21.nml','spine-head-ground-truth-23.nml',...
    'spine-head-ground-truth-24.nml','spine-head-ground-truth-25.nml'});

% load train set
curNodes = table();
gt = struct;
gt.class = categorical({'spinehead'});
gt.segId = [];
gt.label = [];
posSegIds = [];
negSegIds = [];

curMask = gt.class == 'spinehead';
for i = idxTrain
    curNml = slurpNml(nmlFiles{i});
    curNodes = NML.buildNodeTable(curNml);

    curNodes.coord = curNodes.coord + 1;
    curNodes.segId = Seg.Global.getSegIds(param, curNodes.coord);
    assert(all(curNodes.segId));

    curBox = curNml.parameters.userBoundingBox;
    curBox = { ...
        curBox.topLeftX, curBox.topLeftY, curBox.topLeftZ, ...
        curBox.width, curBox.height, curBox.depth};
    curBox = cellfun(@str2double, curBox);

    curBox = Util.convertWebknossosToMatlabBbox(curBox);
    curseg = loadSegDataGlobal(param.seg, curBox);

    posSegIds = reshape(unique(curNodes.segId), [], 1);
    negSegIds = reshape(setdiff(curseg, [0; posSegIds]), [], 1);

    gt.segId = cat(1,gt.segId, double(posSegIds), double(negSegIds));
    labelTemp = zeros(numel(posSegIds)+numel(negSegIds), numel(gt.class));
    labelTemp(1:numel(posSegIds), curMask) = +1;
    labelTemp((numel(posSegIds) + 1):end, curMask) = -1;

    gt.label = cat(1,gt.label, labelTemp);

    clear curNml
end

Util.log('Finding cubeIds...')
cubeIdsSH = cubeIdsForSegs(param.saveFolder, gt.segId(:));
cubeIdsTypeEM = cubeIdsForSegs(param.saveFolder, gtType.segId(:));
cubeIds = unique(vertcat(cubeIdsSH(:), cubeIdsTypeEM(:)));
save(fullfile(nmlDir,'cubeIds.mat'), 'cubeIds')

function cubeIds = cubeIdsForSegs(rootDir, segIds)
    % load meta data
    metaFile = fullfile(rootDir, 'segmentMeta.mat');
    meta = load(metaFile, 'segIds', 'cubeIdx');
    
    [found, row] = ismember(segIds, meta.segIds);
    
    assert(all(found));
    cubeIds = meta.cubeIdx(row);
end
