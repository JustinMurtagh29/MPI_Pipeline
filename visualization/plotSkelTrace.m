function plotSkelTrace( trace, colorString )
%plotSkelTrace(trace, colorString)
toPlot = trace.nodes(trace.edges',1:3);
toPlot = reshape(toPlot,[2,size(trace.edges,1),3]);
plot3(toPlot(:,:,1),toPlot(:,:,2),toPlot(:,:,3), colorString);
end

