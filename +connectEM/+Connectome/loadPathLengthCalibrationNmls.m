function calibT = loadPathLengthCalibrationNmls(param, nmlDir)
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    nmlFiles = dir(fullfile(nmlDir, '*.nml'));
    nmlFiles = reshape({nmlFiles.name}, [], 1);
    nmlFiles = fullfile(nmlDir, nmlFiles);
    
    calibT = table;
    calibT.nmlFile = nmlFiles;
    
    calibT.id = regexp( ...
        calibT.nmlFile, '.*-(?<id>\d+)\.nml$', 'tokens', 'once');
    calibT(cellfun(@isempty, calibT.id), :) = [];
    
    calibT.id = cat(1, calibT.id{:});
    calibT.id = cellfun(@str2double, calibT.id);
    
    % Calculate ground truth path length
    calibT.pathLength = cellfun( ...
        @(nmlFile) calculatePathLength(param, nmlFile), calibT.nmlFile);
    
    % Move file path to the right
    % It's so long that it clutters the display output
    calibT = calibT(:, [2:end, 1]);
end

function pathLength = calculatePathLength(param, nmlFile)
    nml = slurpNml(nmlFile);
    trees = NML.buildTreeTable(nml);
    
    nodes = NML.buildNodeTable(nml);
    nodes.coord = nodes.coord .* param.raw.voxelSize;
    
    % Find tree with oldest node
    % This way we identify the template tracing
   [~, templateTreeId] = min(nodes.time);
    templateTreeId = nodes.treeId(templateTreeId);
    
    calibTreeId = setdiff(nodes.treeId, templateTreeId);
    assert(isscalar(calibTreeId));
    
    calibEdges = trees.edges{trees.id == calibTreeId};
    calibEdges = cat(2, calibEdges.source, calibEdges.target);
   [~, calibEdges] = ismember(calibEdges, nodes.id);
    
    pathLength = ...
        nodes.coord(calibEdges(:, 1), :) ...
      - nodes.coord(calibEdges(:, 2), :);
    pathLength = sqrt(sum(pathLength .^ 2, 2));
    pathLength = sum(pathLength);
end
