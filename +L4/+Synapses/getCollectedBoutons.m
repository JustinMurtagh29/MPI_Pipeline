function [b2agglo, agglo2b] = getCollectedBoutons( agglos, boutons )
%GETCOLLECTEDBOUTONS Get the boutons collected by an axon agglomeration.
% INPUT agglos: [Nx1] cell
%           Cell array of axon agglomerates, i.e. each cell contains an
%           integer list of segment ids belonging to the corresponding axon
%           agglomerate (must be non-overlapping).
%       boutons: [Nx1] cell
%           Cell array of bouton agglomerates.
%           (see also L4.Synapses.boutonAgglo)
% OUTPUT b2agglo: [Nx1] int
%           For each bouton it ontains the linear index of agglos that
%           overlaps with the corresponding bouton. If an bouton overlaps
%           with multiple agglos then it is associated to the agglo with
%           most overlap in terms of number of overlapping segments.
%        agglo2b: [Nx1] cell
%           Cell array of same length as agglos containing the linear
%           indices of the boutons belonging to the corresponding axon
%           agglo.
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

% get a map specifying if an id is in agglos
axIds = cell2mat(agglos);
aggloIdx = repelem((1:length(agglos))', cellfun(@length, agglos));
m = max(cellfun(@max, boutons));
aggloId2Idx = zeros(max([m; axIds]), 1, 'uint32');
aggloId2Idx(axIds) = aggloIdx;

% map bouton segments to agglo idx (keep only mode)
b2agglo = cellfun(@(x)aggloId2Idx(x), boutons, 'uni', 0);
warning('off', 'all'); % mode warning when input is empty
b2agglo = cellfun(@(x)mode(x(x > 0)), b2agglo);
warning('on', 'all')

if nargout > 1
    [comps, aggloIdx] = Util.group2Cell(b2agglo);
    if aggloIdx(1) == 0
        comps = comps(2:end);
        aggloIdx = aggloIdx(2:end);
    end
    agglo2b = cell(length(agglos), 1);
    agglo2b(aggloIdx) = comps;
end

end

