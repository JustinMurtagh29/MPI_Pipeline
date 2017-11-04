function evaluate(chiasmata, agglos)
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>

    % count number of nodes
    chiasmaNodeCounts = cellfun(@(s) sum(s.isIntersection), chiasmata);
    chiasmaCounts = cellfun(@(s) numel(s.ccCenterIdx), chiasmata);
    chiasmaTotalCount = sum(chiasmaCounts);

    % TODO(amotta): Use table instead
    chiasmaSizes = cell2mat(cellfun(@(s) ...
        cellfun(@numel, s.ccNodeIdx(:)), chiasmata, 'Uni', false));
    chiasmaNrExits = cell2mat(cellfun(@(c) ...
        cellfun(@(p) size(p, 1), c.position), chiasmata, 'Uni', false));

    if isfield(agglos, 'solvedChiasma')
        chiasmaSolved = cell2mat(arrayfun(@(c, a) arrayfun( ...
            @(idx) double(a.solvedChiasma(idx)), c{1}.ccCenterIdx(:)), ...
            chiasmata, agglos, 'Uni', false));
        chiasmaSolved = logical(chiasmaSolved);
    else
        chiasmaSolved = false(size(chiasmaNrExits));
    end

    %% report
    fprintf('\nTotal number of chiasmata: %d\n', chiasmaTotalCount);

    temp = sort(chiasmaCounts, 'descend');
    fprintf('\nLargest number of chiasmata per agglomerate:\n');
    fprintf('%s\n', strjoin(arrayfun(@num2str, temp(1:10), 'Uni', false), ', '));

    temp = sort(chiasmaSizes, 'descend');
    fprintf('\nLargest chiasmata (in number of nodes):\n');
    fprintf('%s\n', strjoin(arrayfun(@num2str, temp(1:10), 'Uni', false), ', '));

    temp = sort(chiasmaNodeCounts, 'descend');
    fprintf('\nLargest number of chiasma nodes per agglomerate:\n');
    fprintf('%s\n', strjoin(arrayfun(@num2str, temp(1:10), 'Uni', false), ', '));

    temp = table;
    temp.nrExits = (4:max(chiasmaNrExits))';
    temp.nrChiasmata = arrayfun(@(c) ...
        sum(chiasmaNrExits == c), temp.nrExits);
    temp.nrMarkedSolved = arrayfun(@(c) ...
        sum(chiasmaSolved(chiasmaNrExits == c)), temp.nrExits);
    fprintf('\n'); disp(temp);
end