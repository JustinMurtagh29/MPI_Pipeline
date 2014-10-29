function newCoord = correctNodePositionToBBox( skel, idx, sizeCube )
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here
oldCoord = skel.nodesNumDataAll(idx,3:5);
% Check for problems 
e = find(skel.edges == idx);
[nrEdge, startEnd] = ind2sub(size(skel.edges), e);
if any(size(nrEdge) ~= [1 1]) || any(size(startEnd) ~= [1 1])
    if isempty(nrEdge)
        warning('Disconnected node, removed from skeleton');
        newCoord = [];
        return;
    else
        warning('More than one edge for node outside bbox, removed from skeleton');
        newCoord = [];
        return;
    end
end

for i=1:length(nrEdge)
    if startEnd(i) == 1
        oldCoordBefore = skel.nodesNumDataAll(skel.edges(nrEdge(i),2),3:5);
    else
        oldCoordBefore = skel.nodesNumDataAll(skel.edges(nrEdge(i),1),3:5);
    end
    diff = oldCoord - oldCoordBefore;
    
    whichDimTooSmall = find(oldCoord < 1);
    whichDimTooBig = find(oldCoord > sizeCube);
    if ~isempty(whichDimTooSmall)
        howMuch = abs(oldCoord(whichDimTooSmall)-1)/diff(whichDimTooSmall);
    end
    if ~isempty(whichDimTooBig)
        howMuch = (sizeCube-oldCoord(whichDimTooBig))/diff(whichDimTooBig);
    end
    newCoord(i,:) = oldCoord + howMuch .* diff;
end
newCoord = round(mean(newCoord,1));
end

