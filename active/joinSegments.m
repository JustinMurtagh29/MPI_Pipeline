function [seg, edgesNew, pNew, borderNew, edgesOld, pOld, borderOld] = joinSegments(seg, edges, p, border, lowerCut, upperCut, lowerBounds, upperBounds)

% Structuring element
se = zeros(3,3,3);
se(:,2,2) = 1;
se(2,:,2) = 1;
se(2,2,:) = 1;

% Check for self references
if any(edges(:,1) == edges(:,2))
	error('Self references present in graph');
end

% Logical indices of accepted edges
accepted = p > upperCut;

% Get all accepted edges 
toJoin = edges(accepted,:);

% Find all connected components of objects to join
cc = findCC(toJoin);

% Relabel all objects belonging to a CC to first ID of CC
for i=1:length(cc)
	for j=1:length(cc{i})
		edges(edges == cc{i}(j)) = cc{i}(1);
		seg(imclose(seg == cc{i}(j), se)) = cc{i}(1);
	end
end

% Find self references (edges between objects that are now joined)
idxSelf = edges(:,1) == edges(:,2);
idsSelf = edges(find(idxSelf),1);
uniqueIdsSelf = unique(idsSelf);

% Remove self references
edges(idxSelf,:) = [];
p(idxSelf) = [];
border(idxSelf) = [];

% Join all borders that belong to the same object pair and touch keeping maximum probability
[edgesOld, pOld, borderOld] = joinBorders(uniqueIdsSelf, edges, p, border, size(seg));

% Remove edges with probability below lower cutoff from query (not from possible neighbours)
rejected = pOld < lowerCut;
edgesNew = edgesOld(~rejected,:);
pNew = pOld(~rejected);
borderNew = borderOld(~rejected);

% Remove all borders outside inner bbox from query (same)
lowerBounds = lowerBounds([2 1 3]); % this switch is annoying, any nicer solution, continues in writeKnowledgeDB?
upperBounds = upperBounds([2 1 3]);
outsideBbox = any(bsxfun(@lt,cat(1,borderNew(:).Centroid),lowerBounds),2) | any(bsxfun(@gt,cat(1,borderNew(:).Centroid),upperBounds),2);
edgesNew = edgesNew(~outsideBbox,:);
pNew = pNew(~outsideBbox);
borderNew = borderNew(~outsideBbox);

end

