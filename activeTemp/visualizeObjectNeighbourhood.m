function problemSeg = visualizeObjectNeighbourhood( g, seg, raw, object )

[neighbourhood, weights] = g.findNeighbours(object);
colors = jet(length(neighbourhood)+2);
figure('Units', 'normalized', 'outerposition', [0 0 .4545 0.8992], 'Renderer', 'OpenGL');
% Isosurfaces
subplot(1,2,1);
set(gca, 'ButtonDownFcn', @rotateAroundObject);
obj = seg == object;
issf = isosurface(obj, .1);
l = patch(issf);
set(l, 'FaceColor', colors(end,:), 'EdgeColor', 'none', 'FaceAlpha', 1);
hold on;
k = cell(length(neighbourhood),1);
minimum = min(weights);
maximum = max(weights);
for i=1:length(neighbourhood)
    obj = seg == neighbourhood(i);
    issf = isosurface(obj, .1);
    k{i} = patch(issf);
    set(k{i}, 'FaceColor', colors(i+1,:), 'EdgeColor', 'none', 'FaceAlpha', .1);
end
daspect([25 25 12]);
title(['Neighbourhood of object with segmentation ID: ' num2str(object)], 'FontSize', 18);
grid on;
view(3);
camlight('headlight');
lighting phong;
axis off;
% Strength of weight bar plot
subplot(2,2,2);
hold on;
for i = 1:numel(weights)
    bar(i, weights(i), 'FaceColor', colors(i+1,:));
    label{i} = num2str(neighbourhood(i));
end
set(gca, 'XTick', 1:numel(weights), 'XTickLabel', label);
% Overlay raw with segmentation
subplot(2,2,4);
problemSeg = zeros(size(seg));
for i=1:length(neighbourhood)
	problemSeg(seg == neighbourhood(i)) = i;
end
problemSeg(seg == object) = max(problemSeg(:)) + 1;
for f=1:size(raw,3)
    if sum(reshape(seg(:,:,f), 1, numel(seg(:,:,f))) == object) > 0
        hold off;
        imshow(raw(:,:,f), [60 180]);
        hold on;
        temp = label2rgb(problemSeg(:,:,f), colors(2:end,:), 'w');
        himage = imshow(temp, [60 180]);
        set(himage, 'AlphaData', 0.2 );
        pause(.5);
    end
end

end

