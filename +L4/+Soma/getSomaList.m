function somaCen = getSomaList( filePath )
%GETSOMALIST Get the locations of all somata.
% INPUT filePath: string
%           Path to KAMIN_cells.xlsx file.
% OUTPUT somaCen: [Nx3] double
%           Centroids of non-glia soma cells.
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

somaCen = xlsread(filePath, 'C2:E134');
[~, glia] = xlsread(filePath, 'H2:H134');
isGlia = ~strcmp(glia, 'not glia');

somaCen = somaCen(~isGlia, :);

end

