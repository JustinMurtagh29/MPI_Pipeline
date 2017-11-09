function evaluate(chiasmata, agglos)
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>

    % count number of nodes
    chiasmaNodeCounts = cellfun(@(s) sum(s.isIntersection), chiasmata);
    chiasmaCounts = cellfun(@(s) numel(s.ccCenterIdx), chiasmata);
    chiasmaTotalCount = sum(chiasmaCounts);
    
    chiasmaT = table;
    chiasmaT.size = cell2mat(cellfun( ...
        @(s) cellfun(@numel, s.ccNodeIdx(:)), ...
        chiasmata, 'UniformOutput', false));
    chiasmaT.nrExits = cell2mat(cellfun( ...
        @(c) cellfun(@numel, c.queryIdx(:)), ...
        chiasmata, 'UniformOutput', false));

    if isfield(agglos, 'solvedChiasma')
        chiasmaT.isSolved = cell2mat(arrayfun( ...
            @(c, a) double(arrayfun( ...
                @(idx) a.solvedChiasma(idx), c{1}.ccCenterIdx(:))), ...
            chiasmata, agglos, 'UniformOutput', false));
        chiasmaT.isSolved = logical(chiasmaT.isSolved);
    else
        chiasmaT.isSolved = false(size(chiasmaT, 1), 1);
    end

    %% report
    fprintf('\nTotal number of chiasmata: %d\n', chiasmaTotalCount);

    temp = sort(chiasmaCounts, 'descend');
    fprintf('\nLargest number of chiasmata per agglomerate:\n');
    fprintf('%s\n', strjoin(arrayfun(@num2str, temp(1:10), 'Uni', false), ', '));

    temp = sort(chiasmaT.size, 'descend');
    fprintf('\nLargest chiasmata (in number of nodes):\n');
    fprintf('%s\n', strjoin(arrayfun(@num2str, temp(1:10), 'Uni', false), ', '));

    temp = sort(chiasmaNodeCounts, 'descend');
    fprintf('\nLargest number of chiasma nodes per agglomerate:\n');
    fprintf('%s\n', strjoin(arrayfun(@num2str, temp(1:10), 'Uni', false), ', '));
    
    temp = table;
   [temp.nrExits, ~, tempUniRows] = unique(chiasmaT.nrExits);
    temp.nrChiasmata = accumarray(tempUniRows, 1);
    temp.nrMarkedSolved = accumarray(tempUniRows, chiasmaT.isSolved);
    fprintf('\n'); disp(temp);
end