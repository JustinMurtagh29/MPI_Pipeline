function interfaceClassificationForSeg(p, cubeNo, aggloT, saveFeatures, saveInterfaces)
%INTERFACECLASSIFICATIONFORSEG Classify a local segmentation cube and save
% result to synapse file.
% INPUT pFile: Path to segmentation paramter struct.
%       cubeNo: Integer specifying the linear index of the local segmentation
%           cube in p.local.
%       aggloT: (Optional) Double between 0 and 1 specifying the lower
%           threshold on the GP probability above which segments will get
%           merged.
%           (Default: 1 - no prior merging).
%       saveFeatures: (Optional) Logical specifying whether to save the feature
%           matrix. Features are saved to
%           [p.local(i).saveFolder 'interfaceWeights.mat'].
%           (Default: false)
%       saveInterfaces: (Optional) Logical specifying whether to save the
%           interfaces.
%           Interfaces are saved to
%           [p.local(i).saveFolder 'interfaces.mat'].
%           (Default: false)
%
% see also interfaceClassification
%
% NOTE This function expects to find the interface classifier at
%      p.interfaceClassifier.
%
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

%parse inputs
if ~exist('saveFeatures','var') || isempty(saveFeatures)
    saveFeatures = false;
end
if ~exist('saveInterfaces','var') || isempty(saveInterfaces)
    saveInterfaces = false;
end

%loaded parameter file already in p
pCube = p.local(cubeNo);

%load segmentation
%seg = Seg.Local.getSegSmall(pCube, true);
m = load(pCube.segFile);
seg = uint32(m.seg);
% IDs should be converted to global indices.
seg = Seg.Local.localGlobalIDConversion('LocalToGlobal',pCube,seg); 

%load edges and borders
m = load(pCube.edgeFile);
edges = m.edges;
m = load(pCube.borderFile);
borders = m.borders;
%load interface classifier and feature map
%m = load(p.interfaceClassifier);
m=load('/gaba/u/bstaffle/data/Classifier/20150708T153215/20150708T153215Pred.mat'); % Benedikt's Classifier
classifier = m.classificationStruct;
featureMap = m.options.featureMap;
featureMap.voxelSize = p.raw.voxelSize;
featureMap.names(32) = {'BallAverage'}; % for old features
featureMap.names(33) = {'BallAverage'}; % for old features


%load raw
bboxRaw = pCube.bboxSmall + [-featureMap.border,featureMap.border];
raw = readKnossosRoi(p.raw.root,p.raw.prefix,bboxRaw,'uint8');
if all(~raw)
    error('The raw file read from %s is empty.',p.raw.root);
end
raw = single(raw);
% Normalize raw data to match ex145 on which SynapseClassifier is trained
myStd = std(raw(:));
myMean = mean(raw(:));
raw = (raw./myStd)*22;
raw = raw + 122 - myMean*(22/myStd);

%prior agglomeration
if exist('aggloT','var') && ~isempty(aggloT)
    m = load(pCube.probFile);
    prob = m.prob;
    m = load(pCube.segmentFile);
    segments = m.segments;
    eClasses = Graph.findConnectedComponents(edges(prob > aggloT,:));
    [ seg, edges, borders ] = Seg.Local.applySegEquiv( eClasses, seg, edges, borders, segments );
end

%calculate prediction and save result
[scores, X, interfaces] = interfaceClassification( raw, seg, edges, borders, featureMap, classifier );

%scoresReshaped = NaN(size(edges));
%tmp = [borders(:).Area]> 150;
%scoresReshaped(tmp,:)=reshape(scores,[],2);

m = matfile([p.local(cubeNo).saveFolder 'synapses.mat'], 'Writable', true);
%m = matfile(pCube.synapseFile,'Writable',true);
m.scores = scores;
if saveFeatures
    m = matfile([p.local(cubeNo).saveFolder 'interfaceWeights.mat'], 'Writable', true);
    m.X = X;
    if exist('aggloT','var') && ~isempty(aggloT)
        m.aggloT = aggloT;
        tmp = [borders(:).Area] > 150;
        m.edges = edges(tmp,:);
        borders = borders(tmp);
        m.edgeCentroids = reshape([borders(:).Centroid], 3, length(borders))';
    end
end
if saveInterfaces
    m = matfile([p.local(cubeNo).saveFolder 'interfaces.mat'], 'Writable', true);
    m.interfaces = interfaces;
end

end
