function [ff, errors] = loadFromDirs(p, nmlDirs)
    % Written by
    %   Kevin Boergens <kevin.boergens@brain.mpg.de>
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    if ~iscell(nmlDirs)
        nmlDirs = {nmlDirs};
    end
    
    ff = struct;
   [ff.segIds, ff.neighbours, ff.filenames, ...
    ff.nodes, ff.startNode, ff.comments, errors] = ...
        connectEM.lookupNmlMulti(p, nmlDirs, false);
    
    second = @(x) x(2);
    ff.filenamesShort = cellfun( ...
        @(x) second(strsplit(x, '__')), ff.filenames);
    
    ff = structfun( ...
        @(field) reshape(field, [], 1), ...
        ff, 'UniformOutput', false);
end