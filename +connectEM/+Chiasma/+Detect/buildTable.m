function chiasma = buildTable(chiasmata, agglos)
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    chiasma = table;
    
    %% standard fields
    chiasma.aggloId = repelem( ...
        reshape(1:numel(chiasmata), [], 1), ...
        cellfun(@(c) numel(c.ccCenterIdx), chiasmata));
    chiasma.chiasmaId = cell2mat(cellfun( ...
        @(c) reshape(1:numel(c.ccCenterIdx), [], 1), ...
        chiasmata, 'UniformOutput', false));
    chiasma.nrExits = cell2mat(cellfun( ...
        @(c) cellfun(@numel, c.queryIdx), ...
        chiasmata, 'UniformOutput', false));
    
    %% optional fields, which require `agglos`
    if ~exist('agglos', 'var') || isempty(agglos); return; end
    
    if isfield(agglos, 'solvedChiasma')
        % mark solved chiasmata
        chiasma.isSolved = arrayfun(@(a, c) ...
            agglos(a).solvedChiasma(chiasmata{a}.ccCenterIdx(c)), ...
            chiasma.aggloId, chiasma.chiasmaId);
    end
end