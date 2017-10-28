function [lut,aggloSegIds] = buildLUT(agglos,maxSegId)
    % lut = buildLUT(maxSegId, agglos)
    %   Builds a look-up table `lut` where entry `lut(segId)` contains the
    %   agglomerate ID the segment with ID `segId` belongs to. That is,
    %   
    %   ismember(segId, agglos{aggloId}) == true
    %   if and only if lut(segId) == aggloId
    %
    % Written by
    %  Marcel Beining
    

aggloSegIds = cell2mat(arrayfun(@(x) x.nodes(~isnan(x.nodes(:,4)),4),agglos,'uni',0));
if ~exist('maxSegId','var') || isempty(maxSegId)
    maxSegId = max(aggloSegIds);
end
lut = zeros(maxSegId, 1);
lut(aggloSegIds)  = repelem(1:numel(agglos),arrayfun(@(x) numel(x.nodes(~isnan(x.nodes(:,4)),4)),agglos));  
end