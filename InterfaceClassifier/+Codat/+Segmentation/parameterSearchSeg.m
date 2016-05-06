function result = parameterSearchSeg( pT, cluster, skipMode )
%PARAMETERSEARCHSEG Perform segmentation parameter search.
% This function requires SegEM to be on the path.
% INPUT pT: A parameter struct for the training regions.
%       cluster: (Optional) A cluster object.
%                (Default: getCluster('gpu')).
%                (see also Cluster/getCluster.m)
%       skipMode: (Optional) Logical indicating whether files that are
%           already present should not be calculated again.
%           (Default: false)
% OUTPUT result: Struct containing the spit-merger rates for test and
%           training set as well as the parameters for the respective
%           segmentation.
%
% Author: Manuel Berning <manuel.berning@brain.mpg.de>
% Modified by Benedikt Staffler <benedikt.staffler@brain.mpg.de>

if ~exist('cluster','var') || isempty(cluster)
    cluster = getCluster('cpu');
end
if ~exist('skipMode','var') || isempty(skipMode)
    skipMode = false;
end

%set output folder
param.outputFolder = [pT.saveFolder 'segmentationPS'];
if ~exist(param.outputFolder,'dir')
    mkdir(param.outputFolder);
end
param.figureSubfolder = [param.outputFolder filesep 'figures'];
if ~exist(param.figureSubfolder,'dir')
    mkdir(param.figureSubfolder);
end

% Parameters for morphological reconstruction
param.r = [0 1 2];

% Parameters for segmentation
param.algo(1).fun = func2str(@(seg,pars) Codat.Segmentation.watershedSeg_v1_cortex(seg, pars{:}));
param.algo(1).par = {0.02:0.02:0.9 0:50:200};
param.algo(2).fun = func2str(@(seg,pars) Codat.Segmentation.watershedSeg_v2_cortex(seg, pars{:}));
param.algo(2).par = {[0.01:0.01:0.05 0.1:0.05:0.9 0.95:0.01:0.99] 0:50:200};
param.paramCell = getParamCombinations(param.algo);

% Parameter for evaluation
param.nodeThres = [1, 2];

% load skeletons for split-merger metric
for i = 1:length(pT.local)
    param.skel{i} = parseNml(pT.local(i).trainFileLocal);
    param.totalPathLength(i) = getPathLength(param.skel{i});
end

%morphological reconstruction
fprintf('[%s] - Starting morphological reconstruction:\n', datestr(now));
tic;
idx = 1;
inputCell = [];
for localIdx=1:length(pT.local)
   for radiusIdx=1:length(param.r)
       outputFile = [param.outputFolder filesep 'morphRecon-' num2str(localIdx) '-' num2str(radiusIdx)  '.mat'];
       if skipMode && exist(outputFile,'file')
           continue
       end
       inputCell{idx} = {pT.local(localIdx).class, pT.local(localIdx).bboxSmall, param.r(radiusIdx), pT.cnn.range, outputFile};
       idx = idx + 1;
   end
end
if ~isempty(inputCell)
    functionH = @morphologicalReconstruction;
    job = startJob(cluster, functionH, inputCell);
    fprintf('[%s] - Waiting for job results:\n', datestr(now));
    wait(job, 'finished');
end
toc;
clear job functionH inputCell;

%segmentation
fprintf('[%s] - Starting segmentation:\n', datestr(now));
tic;
idx = 1;
inputCell = [];
for localIdx=1:length(pT.local)
   for radiusIdx=1:length(param.r)
       for segAlg=1:size(param.paramCell,2)
           inputFile = [param.outputFolder filesep 'morphRecon-' num2str(localIdx) '-' num2str(radiusIdx)  '.mat'];
           outputFile = [param.outputFolder filesep 'seg-' num2str(localIdx) '-' num2str(radiusIdx) '-' num2str(segAlg) '.mat'];
           if skipMode && exist(outputFile,'file')
               continue
           end
           %idx = sub2ind([length(pT.local) length(param.r) size(param.paramCell,2)], localIdx, radiusIdx, parVar);
           inputCell{idx} = {param.paramCell{segAlg}, inputFile, outputFile};
           idx = idx + 1;
       end
   end
end
if ~isempty(inputCell)
    functionH = @performSegmentation;
    job = startJob(cluster, functionH, inputCell);
    fprintf('[%s] - Waiting for job results:\n', datestr(now));
    wait(job, 'finished');
end
toc;
clear job functionH inputCell;

%evaluation
fprintf('[%s] - Starting evaluation:\n', datestr(now));
tic
idx = 1;
inputCell = cell(0,1);
resultParams = cell(0,3);
for localIdx=1:length(pT.local)
    for radiusIdx=1:length(param.r)
        for segAlg=1:size(param.paramCell,2)
            segFile = [param.outputFolder filesep 'seg-' num2str(localIdx) '-' num2str(radiusIdx)  '-' num2str(segAlg) '.mat'];
            inputCell{idx} = {segFile, param.skel{localIdx}, param.nodeThres};
            resultParams{idx,1} = radiusIdx;
            resultParams(idx,2) = param.paramCell(segAlg);
            resultParams{idx,3} = [param.outputFolder filesep 'seg-' num2str(localIdx) '-' num2str(radiusIdx)  '-' num2str(segAlg) '.mat'];
            idx = idx + 1;
        end
    end
end
resultParams = cell2table(resultParams,'VariableNames',{'morphR','segParams','filename'});
functionH = @evaluateSegWrapper;
job = startJob(cluster, functionH, inputCell, 1);
fprintf('[%s] - Waiting for job results:\n', datestr(now));
wait(job, 'finished');
%collect outputs
out = fetchOutputs(job);
out = vertcat(out{:});
delete(job);
toc

%reformat output sorted by node thresholds
for i = 1:size(out,2)
    smCurve = cell2mat(out(:,i));
    smCurve = array2table(smCurve,'VariableNames',{'MergeFreeLength','SplitFreeLength', 'AveragePathLengthPerObject','InterErrorDistance'});
    smCurve.Properties.VariableUnits = {'um','um','um','um'};
    %use first cube as training, others as test
    smCurveTrain = smCurve(1:end/length(pT.local),:);
    smCurveTest = smCurve(end/length(pT.local) + 1:end,:);
    result.smCurveTrain{i} = smCurveTrain;
    result.smCurveTest{i} = smCurveTest;
    clear smCurve
end

%save parameters for each result
result.paramsTrain = resultParams(1:end/length(pT.local),:);
result.paramsTest = resultParams(end/length(pT.local) + 1:end,:);
result.paramStruct = param;
end


function morphologicalReconstruction(classParam, bbox, r, range, outputFile)

classification = loadClassData(classParam.root,classParam.prefix,bbox);
bm = Codat.Segmentation.morphRecon(classification, r, range);
m = matfile(outputFile);
m.affReconRecon = bm;

end


function performSegmentation(paramCell, inputFile, outputFile, fullGrow)
% For each training region and radius for morphological reconstruction:
% Perform parameter search as defined in parameters pS

if nargin < 4
    fullGrow = false;
end

m = load(inputFile);
affReconRecon = m.affReconRecon;
% segment and evaluate performance for each set of parameter
fun = str2func(paramCell{1});
segmentation = fun(affReconRecon,paramCell{2});

if fullGrow
    segTemp = imdilate(segmentation, ones(3,3,3));
    borders = segmentation == 0;
    segmentation(borders) = segTemp(borders);
end

m = matfile(outputFile);
m.segmentation = segmentation;
end

function sm = evaluateSegWrapper(segFile, skel, nodeThres)
m = load(segFile);
segmentation = m.segmentation;
result = cell(length(nodeThres),1);
sm = cell(1,length(nodeThres));
pL = getPathLength(skel)./1000; %path lenght in um
for i = 1:length(nodeThres)
    result{i} = evaluateSeg(segmentation, skel, nodeThres(i));
    m = pL./max(result{i}.merge.sum,1);
    s = pL./max(result{i}.split.sum,1);
    ied = 1/(1/m + 1/s);
    sm{i} = [m, s, pL./max(result{i}.general.maxNrObjects,1),ied];
end
end

function inputs = getParamCombinations( algo )
    % Loop over each algorithm
    for i=1:length(algo)
        if ~isempty(algo(i).par)
            numParam = length(algo(i).par);
            numParamVariations = cellfun(@length, algo(i).par);
            % Loop over each parameter
            for j=1:prod(numParamVariations)
                [subs{1:numParam}] = ind2sub(numParamVariations, j);
                inputs{i}{j} = {algo(i).fun cellfun(@(x,s)(x(s)), algo(i).par, subs, 'UniformOutput', false)};
            end
        end
    end
    inputs = [inputs{:}];
end

function eval = evaluateSeg( segmentation, skeletons, nodeThres )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

general = struct('equivMatrix', {}, 'maxNrObjects', {});
nodes = cell(length(size(skeletons,2)),1);
split = struct('vec', {}, 'idx', {}, 'obj', {}, 'sum', {});
merge = struct('vec', {}, 'idx', {}, 'obj', {}, 'sum', {});
for i=1:size(segmentation(1),1)
    for j=1:size(segmentation,2)
        for k=1:size(segmentation,3)
            general(i,j,k).maxNrObjects = single(length(unique(segmentation{i,j,k}(:))));
            general(i,j,k).equivMatrix = zeros(size(skeletons,2), single(max(segmentation{i,j,k}(:))));
            for l=1:size(skeletons,2)
                if size(skeletons{l}.nodes,1) > 0
                    nodes{l} = skeletons{l}.nodes(:,1:3);
                    for m=1:size(nodes{l},1)
                        if segmentation{i,j,k}(nodes{l}(m,1), nodes{l}(m,2), nodes{l}(m,3))
                           general(i,j,k).equivMatrix(l,segmentation{i,j,k}(nodes{l}(m,1), nodes{l}(m,2), nodes{l}(m,3))) = ...
                                general(i,j,k).equivMatrix(l,segmentation{i,j,k}(nodes{l}(m,1), nodes{l}(m,2), nodes{l}(m,3))) + 1;
                        end
                    end
                end
            end
            general(i,j,k).equivMatrixBinary = general(i,j,k).equivMatrix >= nodeThres;
            % Calculate Splits
            split(i,j,k).vec = sum(general(i,j,k).equivMatrixBinary,2);
            split(i,j,k).idx = find(split(i,j,k).vec > 1);
            split(i,j,k).obj = cell(length(split(i,j,k).idx),1);
            for m=1:length(split(i,j,k).idx)
                split(i,j,k).obj{m} = find(general(i,j,k).equivMatrixBinary(split(i,j,k).idx(m),:));
            end
            split(i,j,k).sum = sum(split(i,j,k).vec(split(i,j,k).idx)-1);
            % Calculate Merger
            merge(i,j,k).vec = sum(general(i,j,k).equivMatrixBinary,1);
            merge(i,j,k).idx = find(merge(i,j,k).vec > 1);
            merge(i,j,k).obj = cell(length(merge(i,j,k).idx),1);
            for m=1:length(merge(i,j,k).idx)
                merge(i,j,k).obj{m} = find(general(i,j,k).equivMatrixBinary(:,merge(i,j,k).idx(m)));
            end
            merge(i,j,k).sum = sum(merge(i,j,k).vec(merge(i,j,k).idx)-1);
        end
    end
end

eval.general = general;
eval.nodes = nodes;
eval.split = split;
eval.merge = merge;

end

function totalLength = getPathLength( skeleton )

if isempty(skeleton{1}.parameters)
    voxelsize= [1 1 1];
else
    voxelsize = [str2double(skeleton{1}.parameters.scale.x) ...
        str2double(skeleton{1}.parameters.scale.y) ...
        str2double(skeleton{1}.parameters.scale.z)];
end

totalLength = 0;
for sk=1:length(skeleton)
    if ~isempty(skeleton{sk})
        if size(skeleton{sk}.edges,2) == 1
            % skeleton edges seem to switch dimensionality if only
            % one edge is present
            edgePos = skeleton{sk}.nodes(skeleton{sk}.edges(:,1)',1:3);
            totalLength = totalLength + norm((edgePos(1,:)-edgePos(2,:)).*voxelsize);
        else
            for ed=1:size(skeleton{sk}.edges,1)
                edgePos = skeleton{sk}.nodes(skeleton{sk}.edges(ed,:),1:3);
                totalLength = totalLength + norm((edgePos(1,:)-edgePos(2,:)).*voxelsize);
            end
        end
    end
end

end
