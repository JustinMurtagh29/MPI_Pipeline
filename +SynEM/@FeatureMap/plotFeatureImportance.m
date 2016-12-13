function plotFeatureImportance( obj, imp, plotNo )
%PLOTFEATUREIMPORTANCE Plot feature importance.
% INPUT imp: [Nx1] double where N = obj.numFeatures and contains a positive
%           number measuring the importance of each feature (e.g. from the
%           return from predictorImportance).
%       plotNo: (Optional) [Nx1] int array which plots to produce (see
%           code below).
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

if ~exist('plotNo','var') || isempty(plotNo)
    plotNo = 1:4;
end
makePlot = false(4,1);
makePlot(plotNo) = true;

if isrow(imp)
    imp = imp';
end

obj.setSelectedFeat(true(obj.numFeatures,1));
imp = mat2cell(imp,cellfun(@numel,obj.selectedFeat));
imp = cellfun(@(x,y)reshape(x,size(y)),imp,obj.selectedFeat, ...
    'UniformOutput',false);
featTextClass = cellfun(@(x)strsplit(class(x),'.'), ...
    obj.featTexture, 'UniformOutput', false);
featTextClass = cellfun(@(x)x{end},featTextClass,'UniformOutput',false);
featShapeClass = cellfun(@(x)strsplit(class(x),'.'), ...
    obj.featShape,'UniformOutput',false);
featShapeClass = cellfun(@(x)x{end},featShapeClass,'UniformOutput',false);

%determine features of same class
textGroup = zeros(length(featTextClass),1);
shapeGroup = zeros(length(featShapeClass),1);
groupNamesText = cell(0,1);
groupNamesShape = cell(0,1);
while any(~textGroup)
    idx = strcmp(featTextClass(find(~textGroup,1)),featTextClass);
    textGroup(idx) = max(textGroup) + 1;
    groupNamesText(end+1) = featTextClass(find(idx,1));
end
while any(~shapeGroup)
    idx = strcmp(featShapeClass(find(~shapeGroup,1)),featShapeClass);
    shapeGroup(idx) = max(shapeGroup) + 1;
    groupNamesShape(end+1) = featShapeClass(find(idx,1));
end
shapeGroup = shapeGroup + max(textGroup);
groupNames = cat(2,groupNamesText,groupNamesShape);

% %single feature importance
% figure;
% groups = [textGroup; shapeGroup];
% impTmp = cellfun(@(x)x(:),imp,'UniformOutput',false);
% count = 1;
% for i = 1:max(groups)
%     currImp = cell2mat(impTmp(groups == i));
%     barh(count:count+length(currImp)-1,currImp,'EdgeAlpha',0);
%     hold on
%     count = count + length(currImp);
% end

%total feature importance
totalImp = cellfun(@(x)sum(x(:)),imp);
totalImp = accumarray([textGroup; shapeGroup], totalImp);
if makePlot(1)
    figure;
    barh(flip(totalImp(totalImp > 0)));
    h = gca;
    h.YTick = 1:(sum(totalImp > 0));
    h.YTickLabel = flip(groupNames(totalImp > 0));
    title('Feature class importance')
end

%total feature importance sorted
if makePlot(2)
    figure;
    [~,sortInd] = sort(totalImp, 'descend');
    barh(totalImp(sortInd));
    h = gca;
    h.YTick = 1:length(totalImp);
    h.YTickLabel = groupNames(sortInd);
    title('Feature class importance')
end

%subvolume importance
if makePlot(3)
    if obj.numSubvolumes > 0
        figure;
        subvolImp = cellfun(@(x)sum(sum(x,1),3), ...
            imp(1:length(featTextClass)), 'UniformOutput',false);
        subvolImp = cell2mat(subvolImp);

        %combine subsegments of subvolumes
        subvolImp = mat2cell(subvolImp,ones(1,size(subvolImp,1)), ...
            [1 2.*ones(1,obj.numSubvolumes)]);
        subvolImp = cellfun(@(x)sum(x(:)),subvolImp);
        subvolImpGroup = zeros(max(textGroup),1 + obj.numSubvolumes);
        for i = 1:max(textGroup)
            subvolImpGroup(i,:) = sum(subvolImp(textGroup == i,:),1);
        end
        barh(subvolImpGroup,'stacked');
        h = gca;
        h.YTick = 1:(length(groupNamesText));
        h.YTickLabel = groupNamesText;
        title('Subvolume importance')
        subvolNames = arrayfun(@(x)['Subvol', num2str(x)],...
            1:obj.numSubvolumes,'UniformOutput',false);
        legend('Surface',subvolNames{:});
    end
end

%summary statistic importance
if makePlot(4)
    figure;
    sumStatImp = cellfun(@(x)sum(sum(x,2),3)', ...
        imp(1:length(featTextClass)), 'UniformOutput',false);
    sumStatImp = cell2mat(sumStatImp);
    sumStatImpGroup = zeros(max(textGroup),4 + length(obj.quantiles));
    for i = 1:max(textGroup)
        sumStatImpGroup(i,:) = sum(sumStatImp(textGroup == i,:),1);
    end
    barh(sumStatImpGroup,'stacked');
    h = gca;
    h.YTick = 1:(length(groupNamesText));
    h.YTickLabel = groupNamesText;
    title('Summary statistic importance')
    sumStatNames = arrayfun(@(x)['SumStat', num2str(x)],...
        1:(4 + length(obj.quantiles)),'UniformOutput',false);
    legend(sumStatNames{:});
end

end

