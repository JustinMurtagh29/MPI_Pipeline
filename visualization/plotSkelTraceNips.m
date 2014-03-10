function plotSkelTraceNips( raw, skel, colormap )

for i=1:length(skel)
    toPlot = skel{i}.nodes(skel{i}.edges',1:3);
    toPlot = reshape(toPlot,[2,size(skel{i}.edges,1),3]);
    plot3(toPlot(:,:,1),toPlot(:,:,2),toPlot(:,:,3), 'Color', colormap(i,:), 'LineWidth', 2);
    hold on;
end

axis off;
daspect([25 25 12]);
slices = {[384] [384] [1]};
plotOriginalData(double(raw), slices);


end

