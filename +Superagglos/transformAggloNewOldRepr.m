function aggloOld = transformAggloNewOldRepr(aggloNew)
% transforms the new agglo represenation into the old representation
%
% author marcel.beining@brain.mpg.de

if isfield(aggloNew,'nodes')  % new representation
    empt = arrayfun(@(x) isempty(x.nodes),aggloNew);
    aggloOld = cell(numel(aggloNew),1);
    aggloOld(~empt) = cellfun(@(x) x(~isnan(x(:,4)),4),{aggloNew(~empt).nodes},'uni',0);
    aggloOld = aggloOld(:);
else
    aggloOld = aggloNew(:);
end


