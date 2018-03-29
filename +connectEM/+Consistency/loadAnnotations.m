function synT = loadAnnotations(param, nmlPaths)
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    if ~iscell(nmlPaths)
        nmlPaths = {nmlPaths};
    end
    
    synT = cellfun( ...
        @(nmlPath) forNml(param, nmlPath), ...
        nmlPaths, 'UniformOutput', false);
end

function synT = forNml(param, nmlPath)
    nml = slurpNml(nmlPath);
    nodes = NML.buildNodeTable(nml);
    trees = NML.buildTreeTable(nml);

    % Look up segment IDs
    nodes.coord = nodes.coord + 1;
    nodes.segId = Seg.Global.getSegIds(param, nodes.coord);
    assert(all(nodes.segId));
    
    % Group segment IDs
    agglos = accumarray( ...
        nodes.treeId, nodes.segId, [max(trees.id), 1], ...
        @(segIds) {reshape(segIds, [], 1)}, {zeros(0, 1)});
    
    % Parse tree names
    matches = '^Syn (\d+), (\D+)(?: \d+)?$';
    matches = regexp(trees.name, matches, 'tokens', 'once');
    matches = vertcat(matches{:});
    
    trees.synId = cellfun(@str2double, matches(:, 1));
    trees.isPost = ismember(matches(:, 2), {'Spine head'});
    trees = sortrows(trees, {'synId', 'isPost'});
    
    % Sanity checks
    assert(isequal( ...
        trees.synId(1:2:end), ...
        trees.synId(2:2:end)));
    assert(~any(trees.isPost(1:2:end)));
    assert( all(trees.isPost(2:2:end)));
    
    % Build output
    synT = table;
    synT.id = trees.synId(1:2:end);
    synT.preSegIds = agglos(trees.id(1:2:end));
    synT.postSegIds = agglos(trees.id(2:2:end));
end