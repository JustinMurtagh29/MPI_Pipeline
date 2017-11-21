function sagglos = origFields( sagglos )
%ORIGFIELDS Keep only the 'edges' and 'nodes' field present in the original
%form of superagglos.
% INPUT sagglos: [Nx1] struct
%           Superagglo struct array with any number of fields. All fields
%           except 'nodes' and 'edges are removed.
% OUTPUT sagglos: [Nx1] struct
%           The superagglos struct with all fields but 'nodes' and 'edges'
%           removed.
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

fnames = fieldnames(sagglos);
toKeep = cellfun(@(x)any(strcmp(x, {'nodes', 'edges'})), fnames);
toDel = find(~toKeep);
for i = 1:length(toDel)
    sagglos = rmfield(sagglos, fnames(toDel(i)));
end


end

