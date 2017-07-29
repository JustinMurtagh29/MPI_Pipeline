function interfaces = overlappingSubsegments( interfaces )
%OVERLAPPINGSUBSEGMENTS Revert the exclusive subsegment operations.
% INPUT interfaces: struct
%           see output of SynEM.Svg.calculateInterfaces.
% OUTPUT interfaces: struct
%           The input interface struct with voxels in larger subsegments
%           that also appear in smaller ones removed.
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

for i = length(interfaces.subseg):-1:2
    for j = i-1:-1:1
        for k = 1:size(interfaces.subseg{1}, 1)
            for l = 1:2
            interfaces.subseg{i}{k, l} = sort(cat(1, ...
                interfaces.subseg{i}{k, l}, ...
                interfaces.subseg{j}{k, l}));
            end
        end
    end
end


end

