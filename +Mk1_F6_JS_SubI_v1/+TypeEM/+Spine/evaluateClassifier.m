% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
% Modified by
%   Sahil Loomba <sahil.loomba@brain.mpg.de>

% Use trained classifier and manually check prediction scores on a new bbox

clear;
timeStamp = datestr(now,'yyyymmddTHHMMSS');

if Util.isLocal()
    rootDir = 'E:\u\sahilloo\repos\pipelineHuman';
    addpath(genpath('E:\u\sahilloo\repos\amotta\matlab'))
else
    rootDir = '/tmpscratch/sahilloo/data/Mk1_F6_JS_SubI_v1/pipelineRun_mr2e_wsmrnet';
    addpath(genpath('/gaba/u/sahilloo/repos/amotta/matlab/'))
end
addNoise = false;
noiseDev = 0.01;

vxThrTest= true;
featureSetName = 'segmentAgglomerate';

param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;
param.experimentName = 'Mk1_F6_JS_SubI_v1_mrnet_wsmrnet';
maxSegId = Seg.Global.getMaxSegId(param);
import Mk1_F6_JS_SubI_v1.TypeEM.*

info = Util.runInfo();

segmentMeta = load([param.saveFolder 'segmentMeta.mat'], 'voxelCount', 'point');
vxThr = 100;

% load classifier
m = load('/u/sahilloo/Mk1_F6_JS_SubI_v1/type-em/spineClassifier/20200715T135315.mat');
curClassifier = m.classifier;

% load spinehead training data
nmlDir = fullfile(param.saveFolder, ...
     'tracings', 'box-seeded','spine-head-ground-truth', 'v4');
nmlFiles = fullfile(nmlDir, ...
     {'calibration-box-spines-26.nml','calibration-box-spines-27.nml'});

rng(0);
Util.log('Evaluation on test boxes...')

for idxTest = 1:numel(nmlFiles)
    if addNoise
        experimentName = sprintf('box_%d_addNoise-%.3f_testset', idxTest, noiseDev);
    else
        experimentName = sprintf('box_%d_testset', idxTest);
    end
    
    % load test set
    curNml = slurpNml(nmlFiles{idxTest});
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
    
    if vxThrTest
        % throw away segments that are too small
        vxCount = segmentMeta.voxelCount(negSegIds);
        toDel = vxCount < vxThr;
        negSegIds(toDel) = [];
    end

    gtTest = struct;
    gtTest.class = categorical({'spinehead'});
    gtTest.segId = cat(1, double(posSegIds), double(negSegIds));
    gtTest.label = zeros(numel(gtTest.segId), numel(gtTest.class));
    
    curMask = gtTest.class == 'spinehead';
    gtTest.label(1:numel(posSegIds), curMask) = +1;
    gtTest.label((numel(posSegIds) + 1):end, curMask) = -1;
    
    gtTest = TypeEM.GroundTruth.loadFeatures( ...
        param, curClassifier.featureSetName, gtTest);

    if addNoise
        gtTest.featMat = augmentFeatures(gtTest.featMat, noiseDev);
    end

    gtTest.scores = TypeEM.Classifier.apply(curClassifier, gtTest.featMat);
    % without true labels, conversion from scores to probs is not possible
    gtTest.probs = gtTest.scores;

    % Look at segment scores
    curCount = 100;
    className = 'spinehead';
    skels = Debug.inspectTPs(param, curCount, className, gtTest);
    skels.write(fullfile(param.saveFolder,'typeEM','spine', featureSetName,sprintf('%s-%d-all-labels-%s-%s.nml',timeStamp,idxTest, className, experimentName)));
end

function X = augmentFeatures(X, dev)
%AUGMENTFEATURES Add random noise to features based on the feature mean.
% INPUT X: [NxM] float
%           Feature matrix. Rows correspond to instances and columns to
%           features.
%       dev: (Optional) float
%           Fraction of the feature mean that is used as the standard
%           deviation for random noise for the corresponding feature. i.e.
%           for each column the following is done
%           X(:,i) = X(:,i) + randn(length(X(:,i), 1), 1).*mean(X(:,i))*dev
%           (Default: 0.01)
% OUTPUT X: [NxM] float
%           The input features with added random noise.
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

if ~exist('dev', 'var') || isempty(dev)
    dev = 0.01;
end

m = mean(X, 1);
X = bsxfun(@plus, X, bsxfun(@times, randn(size(X), 'like', X), m.*dev));

end
