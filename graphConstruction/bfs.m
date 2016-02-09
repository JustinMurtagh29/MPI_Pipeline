function components = bfs(idx,adjacency_list)

% Catch obvious case
if length(idx) == 1
   components{1} = idx;
   return;
end

components = {};
k = 1;
nonvisited = 1:length(idx);

while ~isempty(nonvisited)
    % This is reached whenever new component is 'walked'
    % Add first node in nonvisited list to k-th component
    components{k} = idx(nonvisited(1));
    l1 = size(components{k},1);
    % Add all neighbours
    components{k} = cat(1, components{k}, adjacency_list{nonvisited(1)});
    l2 = size(components{k},1);
    % Remove first nodes in nonvisted list
    nonvisited(1)=[];
    % As long as component is, keep walking
    while l1~=l2
        % Collect neighbours of all neighbours
        v = [];
        for j=l1+1:l2
            index = find(idx == components{k}(j));
            v = cat(1, v, adjacency_list{index});
            nonvisited(nonvisited == index) = [];
        end
        % Length of component before agglomeration
        l1 = l2;
        % Add only intersect of neighbours with nonvisited
        v = intersect(v, idx(nonvisited));
        if ~isempty(v)
            components{k} = cat(1, components{k}, v);
        end
        % Update component size
        l2 = size(components{k},1);
    end
    % Increase component counter
    k = k + 1;
end

end
