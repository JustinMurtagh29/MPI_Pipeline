classdef supervoxelGraph
	properties
		supervoxel
		supervoxelContacts
	methods
		function sg = supervoxelGraph(p)
			supervoxel = struct('cubeCoords', {}, 'segmentID', {});
			for i=1:size(p.local,1)
				for j=1:size(p.local,2)
					for k=1:size(p.local,3)
						supervoxel
					end
				end
			end
			supervoxelContacts = sparse
		end

end
