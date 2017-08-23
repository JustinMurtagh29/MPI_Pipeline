function [bIdx, b2Agglo] = getBoutonsOntoAgglo( boutonsPS, postsynAgglos)
%GETBOUTONSONTOAGGLO Get all boutons onto target agglomerates.
% INPUT boutonsPS: [Nx1] cell
%           Postsynaptic segment ids for each bouton as a cell array.
%       postsynAgglos: [Nx1] cell
%           Cell array containing all segment ids of one agglo in each
%           cell. The agglo classes should be postsynaptic target classes.
% OUTPUT bIdx: [Nx1] cell
%           Cell array containing the linear indices of boutons onto the
%           corresponding agglo.
%        b2Agglo: [Nx1] cell
%           Linear indices of the postsynaptic processes for the
%           corresponding bouton.
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

boutonsPS = cellfun(@unique, boutonsPS, 'uni', 0);
idToBouton = accumarray(cell2mat(boutonsPS), ...
    repelem(1:length(boutonsPS), cellfun(@length, boutonsPS))', ...
    [], @(x){x});
bIdx = cellfun(@(x)unique(cell2mat(idToBouton(x))), postsynAgglos, ...
    'uni', 0);
b2Agglo = accumarray(cell2mat(bIdx), ...
    repelem(1:length(bIdx), cellfun(@length, bIdx)), [], @(x){x});
b2Agglo(end+1:length(boutonsPS)) = {[]};

end

