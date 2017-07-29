function [ edgeIdx, edgeDir, label ] = loadInhTestSet( dataFolder )
%LOADINHTESTSET Inhibitory test set loading.
% INPUT dataFolder: string
%           Path to SynEM data folder containing the 'InhibitoryTestSet'
%           folder.
% OUTPUT edgeIdx: [Nx1] int
%           Global indices of the edges w.r.t. 20141007T094904
%           segmentation.
%        edgeDir: [Nx1] int
%           Direction of the corresponding edge:
%           '1': First id in corresponding edge belongs to a test set axon.
%           '2': Second id in corresponding edge belongs to a test set
%                axon.
%        label: [Nx1] int
%           Integer ids that mark all interfaces belonging to the same
%           ground truth synapse.
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

inhTestFolder = fullfile(dataFolder, 'InhibitoryTestSet');

% ground truth data
m = load(fullfile(inhTestFolder, 'GroundTruth.mat'));

% get global edges with direction and label
edgeIdx = cell2mat(m.edgeIdx);
edgeDir = cell2mat(m.edgeDir);
label = m.labelAll;

end

