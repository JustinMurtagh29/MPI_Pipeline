function idx = getFeatureIdx( obj, f_no )
%GETFEATUREIDX Get the column indices for the specified feature.
% INPUT f_no: int
%           Number of the features. Features are enumerated in the order
%           they are listed first in obj.featTexture and then in
%           obj.featShape.
% OUTPUT idx: [1xN] int
%           The column indices of the corresponding feature in the output
%           of obj.calculate.
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

offset = sum(cellfun(@(x)sum(x(:)), obj.selectedFeat(1:(f_no -1))));
idx = (1:sum(obj.selectedFeat{f_no}(:))) + offset;

end

