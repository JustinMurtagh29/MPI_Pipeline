function [h,featureImportanceStat] = featureImportancePlot( imp, featureMap,style )
%FEATUREIMPORTANCEPLOT Extract and plot relevant features for a classifier.
% INPUT imp: The output of the predictorImportance function of a ensemble.
%       featureMap: The feature map on which the classifier was trained.
% OUTPUT h: Handle to the plot
%        featureImportanceStat: Structure containing the summed predictor
%                               importance values for each base feature.
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

featureInfo = featureMap.featureInfo;
indices = imp > 0;
usedFeatures = featureInfo(indices);
indices = find(indices);

%define feature structure containing all base features with derived
%features
s = struct;
sumStatFields = {'total','q1','q2','q3','min','max','mean','var','skew','kur','c','s1','s2','r40','r80','r160'};
for k = 1:length(sumStatFields)
    s.(sumStatFields{k}) = 0;
end
featureImportanceStat = struct;
[names,ia] = unique(featureMap.names);
type = featureMap.type(ia);
for k = 1:length(names)
    featureImportanceStat.(names{k}) = s;
end

%sum all features coming from same base feature
for k = 1:length(usedFeatures)
    str = strsplit(usedFeatures{k},'_');
    ind = strcmp(str{1},featureMap.id);
    featureName = featureMap.names{ind};
    featureImportanceStat.(featureName).total = featureImportanceStat.(featureName).total + imp(indices(k));
    switch featureMap.type{ind}
        case 'texture'
            featureImportanceStat.(featureName).(str{2}) = featureImportanceStat.(featureName).(str{2}) + imp(indices(k));
            featureImportanceStat.(featureName).(str{3}) = featureImportanceStat.(featureName).(str{3}) + imp(indices(k));
            if ~strcmp(str{3},'c')
                featureImportanceStat.(featureName).(str{4}) = featureImportanceStat.(featureName).(str{4}) + imp(indices(k));
            end
    end
end

%gather everything into a single matrix and delete zero columns
featureImpMatrix = cellfun(@(x)cell2mat(struct2cell(x)),struct2cell(featureImportanceStat),'UniformOutput',false);
featureImpMatrix = cell2mat(featureImpMatrix');
ind = any(featureImpMatrix,1);
featureImpMatrix(:,~ind) = [];
names(~ind) = [];
type(~ind) = [];

switch style
    case {'h','horizontal'}
        %total feature importance horizontal plot
        figure;
        bar(featureImpMatrix(1,:),'stack'); %plot the total feature importance
        h = gca;
        set(h,'XTicklabel',names);
        set(h,'FontSize',24);

        position_rectangle = get(h,'Position');
        position_rectangle(2) = 0.38;
        position_rectangle(4) = 0.5;
        set(h,'Position',position_rectangle);

        XTickLabel = get(gca,'XTickLabel');
        set(gca,'XTickLabel',' ');
        hxLabel = get(gca,'XLabel');
        set(hxLabel,'Units','data');
        xLabelPosition = get(hxLabel,'Position');
        y = xLabelPosition(2);
        XTick = get(gca,'XTick');
        y=repmat(y,length(XTick),1);
        fs = get(gca,'fontsize');
        hText = text(XTick, y, XTickLabel,'fontsize',fs);
        set(hText,'Rotation',45,'HorizontalAlignment','right');
        title('Feature Importance')
    case {'v','vertical'}
        %total feature importance
        figure;
        barh(featureImpMatrix(1,end:-1:1));
        h = gca;
        set(h,'FontSize',24);
        h.YTick = 1:length(names);
        set(h,'YTickLabel',names(end:-1:1));
        xlabel('Feature importance value');
        title('Feature importance');
        
        %total feature importance sorted
        figure;
        [~,sortInd] = sort(featureImpMatrix(1,1:end));
        barh(featureImpMatrix(1,sortInd));
        h = gca;
        set(h,'FontSize',24);
        box off
        set(h,'YTickLabel',[]);
        h.YTick = 1:length(names);
        set(h,'YTickLabel',names(sortInd));
        xlabel('Feature importance value');
        title('Feature importance');
        
        %subvolume importance
        figure;
        barh(featureImpMatrix(11:13,strcmp(type,'texture'))','stack');
        h = gca;
        h.YTick = 1:length(names(strcmp(type,'texture')));
        set(h,'YTicklabel',names(strcmp(type,'texture')));
        set(h,'FontSize',24);
        title('Subvolume importance');
        legend('interface surface','first subvolume','second subvolume','location','SouthEast');
        
        %summary statistics importance
        figure;
        barh(featureImpMatrix(2:10,strcmp(type,'texture'))','stack');
        h = gca;
        h.YTick = 1:length(names(strcmp(type,'texture')));
        set(h,'YTicklabel',names(strcmp(type,'texture')));
        set(h,'FontSize',24);
        title('Summary statistics importance');
        legend('0.25 quantile','0.5 quantile','0.75 quantile','min','max','mean','var','skewness','kurtosis','location','SouthEast');
        
        %rinclude importance
        figure;
        barh(featureImpMatrix(14:16,strcmp(type,'texture'))','stack');
        h = gca;
        h.YTick = 1:length(names(strcmp(type,'texture')));
        set(h,'YTicklabel',names(strcmp(type,'texture')));
        set(h,'FontSize',24);
        title('RInclude importance');
        legend('40 nm','80 nm','160 nm','location','SouthEast');
end

