function aggloOld = transformAggloNewOldRepr(aggloNew)
% transforms the new agglo represenation into the old representation
%
% author marcel.beining@brain.mpg.de

if isfield(aggloNew,'nodes')  % new representation
    aggloOld = cellfun(@(x) x(~isnan(x(:,4),4)),{aggloNew.nodes},'uni',0);
    aggloOld = aggloOld(:);
else
    aggloOld = aggloNew(:);
end


