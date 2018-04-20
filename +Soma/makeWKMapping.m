function makeWKMapping( p, outName, outFolder )
%MAKEWKMAPPING Create soma wk mapping.
% INPUT p: struct
%           Segmentation parameter struct.
%       outName: (Optional) string
%           Name of output file.
%       outFolder: (Optional) string
%           Save path.
%           (Default: '.')
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

maxId = Seg.Global.getMaxSegId(p);
m = load(p.agglo.somaFile);
somaAgglos = m.somaAgglos(:, 1);
zeroComp = setdiff(0:maxId, cell2mat(somaAgglos))';
components = cat(1, {zeroComp}, somaAgglos);

if ~exist('outName', 'var') || isempty(outName)
    [~, outName, ~] = fileparts(p.agglo.somaFile);
end
if ~exist('outFolder', 'var') || isempty(outFolder)
    outFolder = '.';
end
Util.log('Saving mapping file to %s.', fullfile(outFolder, outName));
WK.makeWKMapping(components, outName, outFolder);


end
