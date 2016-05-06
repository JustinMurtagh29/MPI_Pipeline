function [ h, target, pred, input, rp, cubeNames ] = ndInterfacePRCurve( folderPath, classificationStruct,varargin )
%NDINTERFACEPRCURVE Plot pr curve for non-direction sensitive interfaces.
% Plot the pr curve for the interfaces by classifying each interface and
% its direction-inverted version and use the logical and result of these
% two classifications for the final classification of this interface.
% INPUT folderPath: Path to the feature folder.
%       classificationStruct: Struct returned by getPredictionVariables
%       varargin{1}: The indices of the learners to use (e.g. 1:10 for the
%                    first 10 weak learners). [] corresponds to using all
%                    learners.
%       varargin{2}: Logical vector indicating which features to use.
    
%parse varargins
if ~isempty(varargin)
    if length(varargin) == 2
        keepFeatures = varargin{2};
    end
    if isempty(varargin{1})
        learners = 1:length(classificationStruct.W);
    else
        learners = varargin{1};
    end
else
    learners = 1:length(classificationStruct.W);
end

%loading and concatenating features and labels
files = what(folderPath);
L = length(files.mat);
input = cell(L,1);
target = cell(L,1);
pred = cell(L,1);
synScores = zeros(0);
y = false(0);
for l = 1:L
    m = matfile([folderPath,filesep,files.mat{l}]);
    X = m.X;
    classLabels = m.classLabels;
    featureMap = m.featureMap;
    
    %area cutoff
    areaIdx = cellfun(@(x)strcmp('area_f1',x),featureMap.featureInfo);
    sizeCutoff = X(:,areaIdx) > 150;
    X = X(sizeCutoff,:);
    classLabels = classLabels(sizeCutoff);
    
    if exist('keepFeatures','var')
        X = X(:,keepFeatures);
    end
    input{l} = X;
    target{l} = classLabels;
    y = cat(1,y,classLabels(1:end/2) | classLabels(end/2 + 1:end));
    [~,scores] = predict(classificationStruct,input{l},learners);
    synScores = cat(1,synScores,scores(classLabels));
    pred{l} = scores;
end
cubeNames = files.mat;

%calculate rp values
synScores = sort(synScores);
rp = zeros(length(synScores),2);
for s = 1:length(synScores)
    yPred = false(0);
    for l = 1:L
        yPred_part = pred{l} >= synScores(s);
        yPred_part = yPred_part(1:end/2) | yPred_part(end/2 + 1:end);
        yPred = cat(1,yPred,yPred_part);
        
    end
    [~,precision,recall] = confusionMatrix(y,yPred);
    rp(s,:) = [recall,precision];
end

%plot pr curve
figure;
scatter(rp(:,1),rp(:,2),'fill');
h = gca;
set(h,'FontSize',24,'FontName','Arial');
xlabel('Recall');
ylabel('Precision');
axis([0,1,0,1]);
grid on
box off

end

