function mask = balanceClasses( y, mask, maxRatio )
%BALANCECLASSES Create target mask to balance training instances.
% INPUT y: Label cube.
%       mask: A priori mask to labels.
%       maxRatio: Maximal ratio betweeneach pair of class (i.e. if
%                 there should be at most 5 times as many points from the
%                 largest class than from the smallest class maxRatio = 5).
% OUTPUT targetMask: Updated mask which is also balanced according to the
%                    label ratio.
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

y = gather(y);
y_unmasked = y(~mask);
labels = unique(y_unmasked);
if length(labels) == 1
    mask = true(size(mask));
else
    nLabels = histcounts(y_unmasked,length(labels));
    [~,idx] = sort(nLabels);
    for i = 2:length(idx)
        ratio = nLabels(idx(i))/nLabels(idx(1));
        if ratio > maxRatio
            numPoints = min(nLabels(idx(1))*maxRatio,nLabels(idx(i)));
            selectedPoints = [false(numPoints,1);true(nLabels(idx(i)) - numPoints,1)];
            selectedPoints = selectedPoints(randperm(length(selectedPoints)));
            points = y == labels(idx(i)) & ~mask;
            points(points) = selectedPoints;
            mask(points) = true;
        end
    end
end
end

