function interfaces = exclusiveSubsegs( interfaces )
%EXCLUSIVESUBSEGS Make subsegments exclusive.
% INPUT interfaces: struct
%           see output of SynEM.Svg.calculateInterfaces.
% OUTPUT interfaces: struct
%           The input interface struct with voxels in larger subsegments
%           that also appear in smaller ones removed.
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

% only keep additional indices for larger rincludes
for i = length(interfaces.subseg):-1:2
    for j = i-1:-1:1
        for k = 1:size(interfaces.subseg{1}, 1)
            for l = 1:2
            interfaces.subseg{i}{k, l} = setdiff(interfaces.subseg{i}{k, l}, ...
                interfaces.subseg{j}{k, l});
            end
        end
    end
end

% store lower rincludes as logical indices for high rincludes - higher
% memory than the solution above in the test example
% for i = 1:(length(interfaces.subseg) - 1)
%     for j = 1:size(interfaces.subseg{i}, 1)
%         for k = 1:2
%             interfaces.subseg{i}{j, k} = ismember( ...
%               interfaces.subseg{i+1}{j, k}, interfaces.subseg{i}{j, k});
%         end
%     end
% end


end

