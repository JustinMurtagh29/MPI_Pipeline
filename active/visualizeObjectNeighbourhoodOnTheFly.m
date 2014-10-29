function visualizeObjectNeighbourhoodOnTheFly( raw, seg, edges, weights, objectID )

colors = distinguishable_colors(length(edges)+1);
figure('position', [1 41 1600 784], 'Renderer', 'OpenGL');
% Isosurfaces
subplot(2,2,1);
obj = seg == objectID;
issf = isosurface(obj, .1);
l = patch(issf);
set(l, 'FaceColor', colors(end,:), 'EdgeColor', 'none', 'FaceAlpha', 1);
hold on;
k = cell(length(edges),1);
for i=1:length(edges)
    obj = seg == edges(i);
    issf = isosurface(obj, .1);
    k{i} = patch(issf);
    set(k{i}, 'FaceColor', colors(i,:), 'EdgeColor', 'none', 'FaceAlpha', .1);
end
daspect([25 25 12]);
title(['Neighbourhood of object with segmentation ID: ' num2str(objectID)], 'FontSize', 18);
grid on;
view(3);
camlight('headlight');
lighting phong;
axis off;
% Strength of weight bar plot
binSize = 0.2;
limits = [-1 1];
x = limits(1)+binSize/2:binSize:limits(2)-binSize/2;
histogram = zeros(length(x),numel(weights));
meanVal = zeros(1,numel(weights));
for i = 1:numel(weights)
    histogram(:,i) = hist(weights{i}, x);
    histogram(:,i) = histogram(:,i) ./ sum(histogram(:,i));
    meanVal(i) = mean(weights{i});
    label{i} = num2str(edges(i));
end
% Histogram plot
subplot(4,2,2);
hold on;
h = bar(x, histogram);
for i=1:size(h,2)
    set(h(i),'facecolor',colors(i,:),'edgecolor','none');
end
xlim([limits(1) limits(2)]);
xlabel('Affinity');
% Mean Value plot
subplot(4,2,4);
hold on;
for i=1:length(meanVal);
    h(i) = bar(i, meanVal(i));
end
for i=1:size(h,2)
    set(h(i),'facecolor',colors(i,:),'edgecolor','none');
end
set(gca, 'XTick', 1:length(meanVal));
set(gca, 'XTickLabel', label);
% Display raw and segmentation
problemSeg = zeros(size(seg));
for i=1:length(edges)
	problemSeg(seg == edges(i)) = i;
end
problemSeg(seg == objectID) = max(problemSeg(:)) + 1;
% Trim problemSeg
obj = seg == objectID;
x = find(any(any(obj,2),3),1,'first'):find(any(any(obj,2),3),1,'last');
y = find(any(any(obj,1),3),1,'first'):find(any(any(obj,1),3),1,'last');
z = find(any(any(obj,1),2),1,'first'):find(any(any(obj,1),2),1,'last');
[x,y,z] = resizeBBox(x,y,z,size(problemSeg), [80 80 20]);
problemSeg = problemSeg(x, y, z);
problemRaw = raw(x,y,z);
for f=1:size(problemSeg,3)
    subplot(6, ceil(size(problemSeg,3)/3), 3*ceil(size(problemSeg,3)/3) + f);
    imshow(problemRaw(:,:,f), [60 180]);
    hold on;
    temp = label2rgb(problemSeg(:,:,f), colors, 'w');
    himage = imshow(temp, [60 180]);
    set(himage, 'AlphaData', 0.3 );
end
set(gcf, 'PaperPositionMode', 'manual', 'PaperType', 'A4', 'PaperUnits', 'normalized', 'PaperOrientation', 'landscape', 'PaperPosition', [0 0 1 1]);
saveas(gcf, ['C:/Users/mberning/Desktop/active/output/segID' num2str(objectID, '%.5i') '.pdf']);
save(['C:/Users/mberning/Desktop/active/output/segID' num2str(objectID, '%.5i') '.mat'], 'problemSeg', 'problemRaw');
close all;
end

