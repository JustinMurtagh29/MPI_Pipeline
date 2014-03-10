function visualizeTasksSupervoxel(seg, v, parameter, i)

colors = jet(length([v(1).missions(i).possibleEnds(:).probability])+1);
border = [250 250 100];
bbox = [v(1).missions(i).errorCenter - border; v(1).missions(i).errorCenter + border];
bbox = bbox - repmat(parameter.bboxBig(:,1) - 1,1,2)';
seg = seg(bbox(1,1):bbox(2,1),bbox(1,2):bbox(2,2),bbox(1,3):bbox(2,3));
% Resort possible ends according to proabability
[~, resortedIdx] = sort([v(1).missions(i).possibleEnds(:).probability], 'descend');
v(1).missions(i).possibleEnds = v(1).missions(i).possibleEnds(resortedIdx);
% Visualization
figure('Units', 'centimeters', 'Position', [0 0 29.7 21], ...
    'PaperType', 'A4', 'PaperUnits', 'centimeters', 'PaperPosition', [0 0 29.7 21], 'PaperOrientation', 'portrait');
hold on;
% Plot start object
startObj = seg == v(1).missions(i).start.id;
startPos = regionprops(startObj, 'Centroid');
startPos = startPos.Centroid;
startObjSmooth = smooth3(startObj, 'gaussian', 5, 1.2);
issf = isosurface(startObjSmooth, .1);
k = patch(issf);
set(k, 'FaceColor', 'b', 'EdgeColor', 'none', 'FaceAlpha', .5);
% Plot other objects
for j=1:length(v(1).missions(i).possibleEnds)
    endObj = seg == v(1).missions(i).possibleEnds(j).id;
    endPos = regionprops(endObj, 'Centroid');
    endPos = endPos.Centroid;
    border = imdilate(startObj, ones(3,3,3))& ~startObj & imdilate(endObj, ones(3,3,3)) & ~endObj;
    borderPos = regionprops(border, 'Centroid');
    borderPos = borderPos.Centroid;
    endObj = smooth3(endObj, 'gaussian', 5, 1.2);
    issf = isosurface(endObj, .1);
    l = patch(issf);
    prob = v(1).missions(i).possibleEnds(j).probability;
    set(l, 'FaceColor', [1-prob prob 0], 'EdgeColor', 'none', 'FaceAlpha', .5, 'Visible', 'off');
    plot3([startPos(1) borderPos(1) endPos(1)], [startPos(2) borderPos(2) endPos(2)], [startPos(3) borderPos(3) endPos(3)], 'Color', [1-prob prob 0], 'MarkerSize', 10, 'LineWidth', 3);
end
view(3);
camlight('headlight');
lighting phong;
daspect([28 28 11.24]);
close all;
end
