function [ X, y, featureMap, cubeNames, intLookupTable ] = ...
    getTrainingDataFrom( folder, type, groupOutput )
%GETTRAININGDATAFROM Load all interface feature files in specified folder.
% Loads all features in a way specified in type. It also applies an area
% thershold of 150 for each interface.
% INPUT folder: string or cell array
%           String or cell array of strings with paths to .mat files. 
%       	 All mat-files in specified folders will be loaded.
%       type: (optional) string
%       	Specify the type of training data.
%             'single': Consider every interface only once.
%             'augment': Consider every interface twice with
%                        direction-inverted subsegments. If an interface is
%                        synaptic, then also its inverted form is synaptic.
%             'direction': (default) Consider every interface twice with
%                          direction-inverted subsegments. An interface
%                          will only be labeled synaptic if the first
%                          subsegments is synaptic.
%       groupOutput: (optional) logical
%           The outputs X and y will be cell arrays containing X and y for
%           each loaded cube.
%           (Default: false)
% OUTPUT X: [NxM] numerical
%           Feature matrix. Each row corresponds to one interface, each
%           column to one feature.
%        y: [Nx1] logical
%           Labels for interfaces in X.
%        featureMap: SynEM.FeatureMap object
%           The featureMap used to calculate X.
%        cubeNames: [Nx1] cell array
%           The names of the files loaded.
%        intLookupTable: Table with file and number of interface in X and y
%           The first row is the index of the matfile in what(folderPath)
%           and the second row the interface index and may thus change if
%           more files are added to the folder.
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

import SynEM.Util.addFilesep
import SynEM.Util.lsFull

if ~exist('type','var') || isempty(type)
    type = 'direction';
end
if ~exist('groupOutput','var') || isempty(groupOutput)
    groupOutput = false;
end

%load and concatenate features and labels
if ischar(folder)
    files = cellstr(lsFull([addFilesep(folder) '*.mat']));
else
    if isrow(folder)
        folder = folder';
    end
    files = cellfun(@(x)cellstr(lsFull([addFilesep(x) '*.mat'])), ...
        folder, 'UniformOutput',false);
    files = vertcat(files{:});
end
L = length(files);
input = cell(L,1);
labels = cell(L,1);
cubeNames = cell(L,1);
intLookup = zeros(0,2,'uint16');
fprintf('[%s] Loading data from %d files.\n',datestr(now), length(files));
for l = 1:L
    m = matfile(files{l});
    X = m.X;
    classLabels = m.classLabels;
    featureMap = m.featureMap;
    
    switch type
        case 'single'
            X = X(1:end/2,:);
            classLabels = classLabels(1:end/2) | ...
                          classLabels(end/2 + 1:end);
            lookup(1:length(classLabels),:) = l;
            lookup(1:length(classLabels),2) = 1:length(classLabels);
        case 'augment'
            tmp = classLabels(1:end/2) | classLabels(end/2 + 1:end);
            classLabels = cat(1,tmp,tmp);
            lookup(1:length(classLabels),:) = l;
            lookup(1:length(classLabels),2) = ...
                repmat((1:(length(classLabels)/2))',2,1);
        case 'direction'
            %standard way of saving them
            lookup(1:length(classLabels),:) = l;
            lookup(1:length(classLabels),2) = ...
                repmat((1:(length(classLabels)/2))',2,1);
        otherwise
            error('Unknown type specified in varargin');
    end
    input{l} = X;
    labels{l} = classLabels;
    [~,cubeNames{l}] = fileparts(files{l});
    intLookup = cat(1,intLookup,lookup);
    clear lookup
    fprintf('.');
end
fprintf('\n');

tmp = cell2mat(labels);
if ~groupOutput
    X = cell2mat(input);
    y = tmp;
else
    X = input;
    y = labels;
end

intLookupTable.files = files;
intLookupTable.interfaceIdx = intLookup;


fprintf('[%s] %d synapses and %d non-synaptic interfaces loaded.\n', ...
        datestr(now),sum(tmp),sum(~tmp));

end

