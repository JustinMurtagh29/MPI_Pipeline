function [ X,y,featureMap,cubeNames,intLookupTable ] = getTrainingDataFrom( folder, type )
%GETTRAININGDATAFROM Load all interface feature files in specified folder.
% Loads all features in a way specified in type. It also applies an area
% thershold of 150 for each interface.
% INPUT folder: Path to .mat files. All mat-files in specified folder
%       	 will be loaded.
%       type: Specify the type of training data. 'direction' is default.
%             'single': Consider every interface only once.
%             'augment': Consider every interface twice with
%                        direction-inverted subsegments. If an interface is
%                        synaptic, then also its inverted form is synaptic.
%             'direction': Consider every interface twice with
%                          direction-inverted subsegments. An interface
%                          will only be labeled synaptic if the first
%                          subsegments is synaptic.
% OUTPUT X: Feature matrix. Each row corresponds to one interface, each
%           column to one feature.
%        y: Labels for interfaces in X.
%        featureMap: The featureMap used to calculate X.
%        cubeNames: The names of the files loaded.
%        intLookupTable: Table with file and number of interface in X and y
%           The first row is the index of the matfile in what(folderPath)
%           and the second row the interface index and may thus change if
%           more files are added to the folder.
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

if ~exist('type','var') || isempty(type)
    type = 'direction';
end

%load and concatenate features and labels
files = what(folder);
L = length(files.mat);
input = zeros(0,'single');
labels = false(0);
cubeNames = cell(L,1);
intLookupTable = zeros(0,2,'uint16');
fprintf('[%s] Loading data from %d files.\n', ...
         datestr(now), length(files.mat));
for l = 1:L
    m = matfile([Util.addFilesep(folder), files.mat{l}]);
    X = m.X;
    classLabels = m.classLabels;
    featureMap = m.featureMap;
    
    %area cutoff
    areaIdx = cellfun(@(x)strcmp('area_f1',x),featureMap.featureInfo);
    sizeCutoff = X(:,areaIdx) > 150;
    X = X(sizeCutoff,:);
    classLabels = classLabels(sizeCutoff);
    
    switch type
        case 'single'
            X = X(1:end/2,:);
            classLabels = classLabels(1:end/2) | classLabels(end/2 + 1:end);
            lookup(1:length(classLabels),:) = l;
            lookup(1:length(classLabels),2) = 1:length(classLabels);
        case 'augment'
            tmp = classLabels(1:end/2) | classLabels(end/2 + 1:end);
            classLabels = cat(1,tmp,tmp);
            lookup(1:length(classLabels),:) = l;
            lookup(1:length(classLabels),2) = repmat((1:(length(classLabels)/2))',2,1);
        case 'direction'
            %standard way of saving them
            lookup(1:length(classLabels),:) = l;
            lookup(1:length(classLabels),2) = repmat((1:(length(classLabels)/2))',2,1);
        otherwise
            error('Unknown type specified in varargin');
    end
    input = cat(1,input,X);
    labels = cat(1,labels,classLabels);
    cubeNames{l} = files.mat{l};
    intLookupTable = cat(1,intLookupTable,lookup);
    clear lookup
end

X = input;
y = labels;

fprintf('[%s] %d synapses and %d non-synaptic interfaces loaded.\n', ...
        datestr(now),sum(y),sum(~y));

end

