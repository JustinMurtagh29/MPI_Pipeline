function [p, job] = pipelineRun( p, cluster )
%PIPELINERUN Run SynEM as part of the pipeline.
% INPUT p: struct
%           Segmentation parameter struct.
%       cluster: (Optional) parallel.cluster object
%           The cluster object for prediction.
%           (Default: Cluster.getCluster())
% OUTPUT p: struct
%           Segmentation parameter struct updated with the fields
%           p.synEM: Path of classifier and feature map in the
%               segmentation main folder.
%           p.local.synapseFile: Path to the synapse scores in each local
%               segmentation cube folder.
%           p.local.interfaceFeatureFile: Path to the interface features in
%               each local segmentation cube folder.
%        job: parallel.job object
%           The synapse detection cluster job object.
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

%path to classifier relative to this file
curFile = mfilename('fullpath');
[curFolder,~,~] = fileparts(curFile);
curFolder = strsplit(curFolder, filesep);
classifierPath = [fullfile(curFolder{1:end-1}) filesep 'data' filesep ...
    'SynEMPaperClassifier.mat'];
if ~ispc
    classifierPath = [filesep classifierPath];
end

%load classifier
m = load(classifierPath, 'classifier');
classifier = m.classifier;
fm = classifier.options.fm;

%get classifier
if ~exist('cluster', 'var') || isempty(cluster)
    cluster = Cluster.getCluster('-pe openmp 1', ...
                   '-l h_vmem=20G,h_rt=100:00:00', ...
                   '-p -500');
end

%run prediction
[p, job] = SynEM.Seg.predictDataset(p, fm, classifier, [], cluster, ...
    [], true, false);

end

