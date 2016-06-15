function buildAllQueries(param, outDir, edges)
    % buildAllQueries(param, outDir, edges)
    %   Builds all queries from a list of edges and saves
    %   them as NML files to 'outDir'.
    %
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    edgeCount = size(edges, 1);
    
    % loading centre of mass
    m = load([param.saveFolder, 'globalCoMList.mat']);
    com = m.globalCoMList;
    
    % run
    doFunc = @(idx) ...
        writeQuery(outDir, com, edges(idx, :));
    arrayfun(doFunc, 1:edgeCount);
end

function writeQuery(outDir, com, edge)
    fileName = [ ...
        'edge-', num2str(edge(1)), ...
        '-to-', num2str(edge(2)), '.nml'];
    filePath = fullfile(outDir, fileName);
    
    skel = buildEdgeQuery(com, edge);
    skel.write(filePath);
end
