function soma = getSoma(wcAgglos, somaAgglos)
%GETSOMA Get all somatic nodes of whole cell agglos.
%
% INPUT wcAgglos: [Nx1] struct
%           Whole cell agglomeratios in the superagglo format.
%       somaAgglos: [Nx1] cell or struct
%           Soma agglomerates in the agglo or superagglo format.
%
% OUTPUT soma: [Nx1] cell
%           Cell array with the logical indices of somatic nodes for each
%           whole cell agglomerate.
%
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

if iscell(somaAgglos)
   somaIds = cell2mat(somaAgglos);
else
    somaIds = cell2mat(arrayfun(@(x)x.nodes(:,4), somaAgglos, 'uni', 0));
end

soma = cellfun(@(x)ismember(x(:, 4), somaIds), {wcAgglos.nodes}, 'uni', 0);
soma = soma(:);

end