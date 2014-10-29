function [edges, p, border] = joinBorders(uniqueIdsSelf, edges, p, border, bboxSize)
% Handle borders that might have to be joined due to joined objects

for i=1:length(uniqueIdsSelf)
	% Find all unique object ids connected to a joined object
	found = edges == uniqueIdsSelf(i);
	found = found(:, [2 1]);
	found(sum(found,2) == 2,:) = 0;	
	idsOther = edges(found);
	uniqueIdsOther = unique(idsOther);
	for j=1:length(uniqueIdsOther)
		% Determine whether objects have more than one area of touch
		relevantEdges = all(bsxfun(@eq, edges, [uniqueIdsSelf(i) uniqueIdsOther(j)]),2) | ...
				all(bsxfun(@eq, edges, [uniqueIdsOther(j) uniqueIdsSelf(i)]),2);
		idx = find(relevantEdges);
		if length(idx) > 1;
			% Put all borders between pair of objects back into array & check for touch
			temp = zeros(bboxSize, 'uint8');
			for k=1:length(idx)
				temp(border(idx(k)).PixelIdxList) = 1;
			end
			newBorders = regionprops(temp, {'Area' 'Centroid' 'PixelIdxList'});
			% If more than one new border find correspondence of old (border(idx)) and new (newBorders) arrays
			if length(newBorders) > 1
				for k=1:length(newBorders)
					for l=1:length(idx)
						newBordersIdx{k}(l) = any(ismember(border(idx(l)).PixelIdxList, newBorders(k).PixelIdxList));
					end
				end
			else
				newBordersIdx = {1:length(idx)};
			end
			% Add consolidated borders/p/edges
			for k=1:length(newBordersIdx)
				border(end+1) = newBorders(k);
				p(end+1) = max(p(idx(newBordersIdx{k})));
				edges(end+1,:) = [uniqueIdsSelf(i) uniqueIdsOther(j)];
			end
			% Remove old ones
			edges(idx,:) = [];
			p(idx) = [];
			border(idx) = [];	
		end
	end
end	

end

