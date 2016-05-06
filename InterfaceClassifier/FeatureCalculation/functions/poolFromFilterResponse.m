function features = poolFromFilterResponse( filterResponse, interfaceSurfaceList, subsegmentsList, featureIndices, quantiles, summaryStatistics, border )
%POOLFROMFILTERRESPONSE Calculate summary statistics over volumes.
% Calculates the summary statistics specified in featureIndices pooled from
% the filter responseover the interface surfaces and the subsegments.
% Quantiles and other summary statistics can be passed via varargin or
% standard pooling quantities are used which comprise the 0.25,0.5 and 0.75
% quantile as well as the min, max, mean, var, skewness and kurtosis.
% Note: min and max could also be implemented as 0 and 1 quantile which
% seems to be slower in general.
% OUTPUT features: array of size
%                  length(interfaceSurfaceList)x(length(quantiles) +
%                  length(summaryStatistics)), where rows correspond to a
%                  single interface and columns correspond to different
%                  pooling statistics
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

fprintf('%s INFO - pooling from filter response\n',datestr(now,'yyyy-mm-dd HH:MM:SS,FFF'));

%crop border
filterResponse = filterResponse(border(1) + 1:end - border(1), border(2) + 1: end - border(2), border(3) + 1: end - border(3));
summaryStatistics = cellfun(@str2func,summaryStatistics,'UniformOutput',false);
if all(featureIndices(:))
    %this way of calculating all features was quicker in my tests
    interfaceSurfaceFeatures = cellfun(@(x)filterResponse(x),interfaceSurfaceList,'UniformOutput',false);
    singleInterfaceSurfaceFeatures = cellfun(@(x)[quantile(x,quantiles),min(x),max(x),mean(x),var(x),skewness(x),kurtosis(x)],interfaceSurfaceFeatures,'UniformOutput',false);
    tmp = cat(1,singleInterfaceSurfaceFeatures,singleInterfaceSurfaceFeatures);
    features = cell2mat(tmp);
    subsegmentFeatures = cellfun(@(x)filterResponse(x),subsegmentsList,'UniformOutput',false);
    singleSubsegmentFeatures = cellfun(@(x)[quantile(x,quantiles),min(x),max(x),mean(x),var(x),skewness(x),kurtosis(x)],subsegmentFeatures,'UniformOutput',false);
    tmp = cat(1,singleSubsegmentFeatures,singleSubsegmentFeatures(:,reshape([2:2:size(subsegmentsList,2);1:2:size(subsegmentsList,2)],1,size(subsegmentsList,2))));
    tmp = cell2mat(tmp);
    features = cat(2,features,tmp);
else
    %get relevant summary statistics
    quantileIndices = featureIndices(1:length(quantiles),:);
    sumStatIndices = featureIndices(length(quantiles) + 1:end,:);
    
    %preallocate arrays
    interfaceSurfaceFeatures = zeros(length(interfaceSurfaceList),sum(featureIndices(:,1)),'single');
    subsegmentFeatures = cell(1,size(subsegmentsList,2));

    %distribute summary statistics to interface and subsegments
    interfaceQuantiles = quantiles(quantileIndices(:,1));
    interfaceSumStats = summaryStatistics(sumStatIndices(:,1));
    subsegQuantiles = cell(size(subsegmentsList,2),1);
    subsegSumStats = cell(size(subsegmentsList,2),1);
    for iSubseg = 1:size(subsegmentsList,2)
        subsegQuantiles{iSubseg} = quantiles(quantileIndices(:,1 + iSubseg));
        subsegSumStats{iSubseg} = summaryStatistics(sumStatIndices(:,1 + iSubseg));
    end
    
    %pool from filter response
    for k = 1:length(interfaceSurfaceList)
        contactSurfaceResponse = filterResponse(interfaceSurfaceList{k});
        %pool quantiles
        interfaceFeatureRow = [];
        if ~isempty(interfaceQuantiles)
            interfaceFeatureRow = quantile(contactSurfaceResponse,interfaceQuantiles);
        end
        %pool other sumStats
        for l = 1:length(interfaceSumStats)
            interfaceFeatureRow = cat(2,interfaceFeatureRow,interfaceSumStats{l}(contactSurfaceResponse));
        end
        if ~isempty(interfaceFeatureRow)
            interfaceSurfaceFeatures(k,:) = interfaceFeatureRow;
        end

        %features for subsegments
        for iSubseg = 1:size(subsegmentsList,2);
            subsegResponse = filterResponse(subsegmentsList{k,iSubseg});
            subsegFeatureRow = [];
            if ~isempty(subsegQuantiles{iSubseg})
                subsegFeatureRow = cat(2,subsegFeatureRow,quantile(subsegResponse,subsegQuantiles{iSubseg}));
            end
            for l = 1:sum(sumStatIndices(:,1 + iSubseg))
                subsegFeatureRow = cat(2,subsegFeatureRow,subsegSumStats{iSubseg}{l}(subsegResponse));
            end
            subsegmentFeatures{iSubseg} = cat(1,subsegmentFeatures{iSubseg},single(subsegFeatureRow));
        end
    end
    %cat features and duplicate with interchanged subsegments
    features = cat(2,interfaceSurfaceFeatures,cell2mat(subsegmentFeatures(~cellfun(@(x)isempty(x),subsegmentFeatures))));
    subsegmentFeatures = subsegmentFeatures(reshape([2:2:size(subsegmentsList,2);1:2:size(subsegmentsList,2)],1,size(subsegmentsList,2)));
    features = cat(1,features,cat(2,interfaceSurfaceFeatures,cell2mat(subsegmentFeatures(~cellfun(@(x)isempty(x),subsegmentFeatures)))));
end

end

 %quicker calculation of skewness and kurtosis dropping all code checking
 function s = skewness(x)
 x0 = x - mean(x);
 s2 = mean(x0.^2);
 m3 = mean(x0.^3);
 s = single(m3./s2.^(1.5));
 end
 
 function k = kurtosis(x)
 x0 = x - mean(x);
 s2 = mean(x0.^2);
 m4 = mean(x0.^4);
 k = single(m4./s2.^2);
 end
 
 function [ Y ] = quantile( X,p )
 X = sort(X);
 L = length(X);
 Y = zeros(length(p),1);
 for i = 1:length(p)
     a = L*p(i) + 0.5;
     b = floor(a);
     nonint = a - b;
     if nonint == 0
         Y(i) = X(a);
     else
         Y(i) = X(b) + nonint*(X(b + 1) - X(b));
     end
 end
 Y = single(Y');
 end

