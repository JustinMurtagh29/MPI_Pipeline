function [exits, valid] = parseQuery(nmlPath)
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    
    try
        nml = slurpNml(nmlPath);
        nodes = NML.buildNodeTable(nml);
        trees = NML.buildTreeTable(nml);
        comments = NML.buildCommentTable(nml);
    catch
        warning('Invalid NML file %s.', nmlPath);
        
        nodes = table;
        nodes.id = zeros(0, 1);
        nodes.treeId = zeros(0, 1);
        
        trees = table;
        trees.name = cell(0, 1);
        
        comments = table;
        comments.node = zeros(0, 1);
        comments.comment = cell(0, 1);
    end
    
    valid = struct;
    
    % check if there are invalid trees
    invalidTreeName = regexp(trees.name, '^Branch \d+$');
    invalidTreeName = cellfun(@isempty, invalidTreeName);
    valid.nrInvalidTrees = sum(invalidTreeName);
    
    % find tree for each comment
   [~, comments.treeId] = ismember(comments.node, nodes.id);
    comments.treeId = nodes.treeId(comments.treeId);
    
    % parse exit ID
    comments.exitId = regexp( ...
        comments.comment, '^Exit (\d+)', 'tokens', 'once');
    comments.isExit = ~cellfun(@isempty, comments.exitId);
    valid.nrInvalidComments = sum(~comments.isExit);
    
    % split off exits
    exits = table;
    exits.id = comments.exitId(comments.isExit);
    exits.groupId = comments.treeId(comments.isExit);
    
    % fix `exits`
    exits.id = cellfun(@(id) str2double(id{1}), exits.id);
   [~, ~, exits.groupId] = unique(exits.groupId);
    exits = sortrows(exits, 'id');
end