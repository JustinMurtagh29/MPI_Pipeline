function sagglo = toStructOfArray( sagglo )
%TOSTRUCTOFARRAY Convert a superagglo in the array of struct format to a
% struct of arrays format.
% INPUT sagglo: [Nx1] struct
%           Superagglo struct with the fields 'nodes' and 'edges'.
% OUTPUT sagglos: struct
%           Struct with the fields
%           'nodes': [Nx1] cell
%           'edges': [Nx1] cell
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

tmp.nodes = {sagglo.nodes}';
tmp.edges = {sagglo.edges}';
sagglo = tmp;

end

