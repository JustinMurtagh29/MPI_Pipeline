function fp = getFlightPathSegIds2( p, fp )
%GETFLIGHTPATHSEGIDS2 Simple flight path segment id loading.
% INPUT p: struct
%           Segmentation parameter struct.
%       fp: struct
%           Struct containing a field nodes with 3d flight path nodes.
%           (see also Superagglos.getFlightPath)
% OUTPUT fp: struct
%           The input struct with an additional field segIds containing the
%           segment ids for the corresponding nodes in fp.nodes.
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

nodes = {fp.nodes}';
l = cellfun(@(x)size(x, 1), nodes);
if strcmp(p.seg.backend, 'wkwrap')
    segIds = Seg.Global.getSegIds( p, cell2mat(nodes), [1024, 1024, 1024]);
else
    segIds = Seg.Global.getSegIds( p, cell2mat(nodes));
end
segIds = mat2cell(segIds, l, 1);
for i = 1:length(fp)
    fp(i).segIds = segIds{i};
end

end

