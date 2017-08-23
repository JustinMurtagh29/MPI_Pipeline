function p = updateParamsToNewestFiles( p )
%UPDATEPARAMSTONEWESTFILES Update the paths in p to the newest pipeline
%result files.
% INPUT p: struct
%           Segmentation parameter struct.
% OUTPUT p: struct
%           Update segmentation parameter struct.
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

if ~isfield(p, 'svg')
    p = Seg.Util.addSvgFiles(p);
end

graphFile = 'graphNew.mat';
corrFile = 'correspondencesNew2.mat';

p.svg.graphFile = fullfile(p.saveFolder, graphFile);
p.svg.correspondenceFile = fullfile(p.saveFolder, corrFile);

end

