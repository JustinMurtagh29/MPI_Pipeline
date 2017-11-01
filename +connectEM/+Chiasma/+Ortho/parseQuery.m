function [exits, otherComments] = parseQuery(nmlPath)
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    nml = slurpNml(nmlPath);
    nodes = NML.buildNodeTable(nml);
    comments = NML.buildCommentTable(nml);
    
    % find tree for each comment
   [~, comments.treeId] = ismember(comments.node, nodes.id);
    comments.treeId = nodes.treeId(comments.treeId);
    
    % parse exit ID
    comments.exitId = regexp( ...
        comments.comment, '^Exit (\d+)', 'tokens', 'once');
    comments.isExit = ~cellfun(@isempty, comments.exitId);
    
    % split off exits
    exits = table;
    exits.id = comments.exitId(comments.isExit);
    exits.groupId = comments.treeId(comments.isExit);
    
    % fix `exits`
    exits.id = cellfun(@(id) str2double(id{1}), exits.id);
   [~, ~, exits.groupId] = unique(exits.groupId);
    exits = sortrows(exits, 'id');
   
    % find other comments
    otherComments = comments.comment(~comments.isExit);
end