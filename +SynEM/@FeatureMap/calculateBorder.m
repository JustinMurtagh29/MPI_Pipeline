function border = calculateBorder(obj)
%CALCULATEBORDER Calculate the border of the feature map.
% OUTPUT border: [1x3] int containing the calculated obj.border property.
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

obj.border = max(cell2mat( ...
    cellfun(@(x)x.border, obj.featTexture, ...
    'UniformOutput', false)'),[],1);
border = obj.border;
end