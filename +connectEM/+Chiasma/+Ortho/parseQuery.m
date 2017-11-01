function [exits, otherComments] = parseQuery(nmlPath)
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    
    try
        nml = slurpNml(nmlPath);
        nodes = NML.buildNodeTable(nml);
        comments = NML.buildCommentTable(nml);
    catch e
       warning('Invalid NML file %s.', nmlPath);
       
       nodes = table;
       nodes.id = zeros(0, 1);
       nodes.treeId = zeros(0, 1);
       
       comments = table;
       comments.node = zeros(0, 1);
       comments.comment = cell(0, 1);
    end
    
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