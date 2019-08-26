% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
outDir = '/home/amotta/Desktop';

% Configurations
confs = struct;

confs(1).pos = [2910, 4877, 2363] + 1;
confs(1).edge = [10032399, 10032148];

% For more candidate locations, see
% 2019-04-12-ConnectEM-Illustration-for-Supplements/notes.md

fovNm = 2000;
rawRange = [60, 180];

% SynEM's feature map
% Values extracted from +SynEM/data/SynEMPaperClassifier.mat
% Fields of classifier.options.fm 
fmFile = fullfile('+SynEM', 'data', 'SynEMPaperClassifier.mat');

info = Util.runInfo();
Util.showRunInfo(info);

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

fm = load(fullfile(info.git_repos{1}.local, fmFile));
fm = fm.classifier.options.fm;

%% Build figures
for curConfIdx = 1:numel(confs)
    curConf = confs(curConfIdx);
    curPos = curConf.pos;
    curEdge = sort(curConf.edge);
    
    curOutFile = sprintf('connect-em_edge-%d_', curConfIdx);
    curOutFile = fullfile(outDir, curOutFile);
    
    curParam = curPos(:) - param.bbox(:, 1);
    curParam = 1 + floor(curParam ./ param.tileSize)';
    curParam = param.local(curParam(1), curParam(2), curParam(3));
    
    curFovBox = fovNm / param.raw.voxelSize(1);
    curFovBox = [-1, +1] / 2 .* [curFovBox; curFovBox; 0];
    curFovBox = round(curPos(:) + curFovBox);
    
    curBox = curParam.bboxSmall;
    curBoxRel = 1 + curFovBox - curBox(:, 1);
    
    curCrop = @(data) data( ...
        curBoxRel(1, 1):curBoxRel(1, 2), ...
        curBoxRel(2, 1):curBoxRel(2, 2), ...
        curBoxRel(3, 1):curBoxRel(3, 2));
    
    curEdges = load(curParam.edgeFile);
    curEdges = curEdges.edges;
    
    curBorders = load(curParam.borderFile);
    curBorders = curBorders.borders;
    
   	curMask = ismember(curEdges, curEdge, 'rows');
    curEdges = curEdges(curMask, :);
    curBorders = curBorders(curMask);
    
    % Sanity check
    if size(curEdges, 1) ~= 1; continue; end
    if ~isscalar(curBorders); continue; end
    
    %% Build subsegments
    curSeg = loadSegDataGlobal(param.seg, curBox);
    
   [curInterf, curMask] = ...
        SynEM.Svg.calculateInterfaces( ...
            curSeg, curEdges, curBorders, ...
            fm.areaT, param.raw.voxelSize, fm.subvolsSize);
    if ~all(curMask); continue; end
    
    data = struct;
    data.segments = curSeg;
    data.subsegmentsList = curInterf.subseg;
    data.interfaceSurfaceList = curInterf.surface;
    
    %% Build subsegment masks
    % Adapted from +Paper/+SynEM/+Figures/figure1.m
    % of https://gitlab.mpcdf.mpg.de/connectomics/Benedikt
    % Commit 0058c0cf42521f514828840ef26cc20691699fe3
    clmap40 = [8,69,148;0,90,50;255,0,0]/255;
    clmap80 = [49,130,189;49,163,84;255, 0, 0]/255;
    clmap160 = [107,174,214;116,196,118;255, 0, 0]/255;
    L1 = zeros(size(curSeg));
    L2 = zeros(size(curSeg));
    L3 = zeros(size(curSeg));
    L1(data.subsegmentsList{1}{1,1}) = 1;
    L1(data.subsegmentsList{1}{1,2}) = 2;
    L2(data.subsegmentsList{2}{1,1}) = 1;
    L2(data.subsegmentsList{2}{1,2}) = 2;
    L3(data.subsegmentsList{3}{1,1}) = 1;
    L3(data.subsegmentsList{3}{1,2}) = 2;
    L3(data.interfaceSurfaceList{1}) = 3;
    %get rid of overlaps for illustration
    L2(L1 > 0) = 0;
    L3(L1 > 0 | L2 > 0 ) = 0;
    l1 = curCrop(L1);
    l1 = permute(l1,[2 1]);
    l2 = curCrop(L2);
    l2 = permute(l2,[2 1]);
    l3 = curCrop(L3);
    l3 = permute(l3,[2 1]);
    im1 = label2rgb(l1,clmap40);
    im2 = label2rgb(l2,clmap80);
    im3 = label2rgb(l3,clmap160);
    
    %% Load EM and CNN data    
    curRaw = double(loadRawData(param.raw, curFovBox));
    curRaw = (curRaw - rawRange(1)) / (rawRange(2) - rawRange(1));
    
    curClass = loadClassData(param.class, curFovBox);
    curClass = ((curClass / 1.7159) + 1) / 2;
    
    curSeg = loadSegDataGlobal(param.seg, curFovBox);
    
    %% Show figures
    % Adapted from +Paper/+SynEM/+Figures/figure1.m
    % of https://gitlab.mpcdf.mpg.de/connectomics/Benedikt
    % Commit 0058c0cf42521f514828840ef26cc20691699fe3
    curFix = @(d) repmat(permute(d, [2, 1]), [1, 1, 3]);
    
    curFig = figure;
    curAx = axes(curFig); %#ok
    hold(curAx, 'on');
    
    imshow(curFix(curRaw), 'Parent', curAx);
    export_fig(strcat(curOutFile, 'raw.png'), curFig);
    
    
    curFig = figure;    
    curAx = axes(curFig); %#ok
    hold(curAx, 'on');
    
    imshow(curFix(curClass), 'Parent', curAx);
    export_fig(strcat(curOutFile, 'cnn.png'), curFig);
    
    
    curFig = figure;    
    curAx = axes(curFig); %#ok
    hold(curAx, 'on');
    
    imshow(curFix(curClass), 'Parent', curAx);
    himage1 = imshow(im1, 'Parent', curAx);
    set(himage1,'Alphadata',0.7.*(l1 > 0));
    himage2 = image(im2, 'Parent', curAx);
    set(himage2,'Alphadata',0.7.*(l2 > 0));
    himage3 = image(im3, 'Parent', curAx);
    set(himage3,'Alphadata',0.7.*(l3  > 0));
    
    export_fig(strcat(curOutFile, 'subvol.png'), curFig);
    
    
    curFig = figure;
    curAx = axes(curFig); %#ok
    hold(curAx, 'on');
    
    curSize = size(curSeg);
   [curUniSegIds, ~, curUniSeg] = unique(curSeg);
    assert(curUniSegIds(1) == 0);
    
    curColors = distinguishable_colors(numel(curUniSegIds) - 1, [0, 0, 0]);
    curColors = vertcat([0, 0, 0], curColors); %#ok
    
    curSeg = curColors(curUniSeg, :);
    curSeg = reshape(curSeg, [curSize, 3]);
    
    imshow(permute(curSeg, [2, 1, 3]), 'Parent', curAx);
    export_fig(strcat(curOutFile, 'seg.png'), curFig);
end
