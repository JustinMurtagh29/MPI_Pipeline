function plotSkelTraceMinicube( trace, posM, sizeM, colorString )
%plotSkelTraceMinicube(trace, posM, sizeM, colorString)
isWithin = trace.nodes(:,1:3) <= repmat(posM' + sizeM, [size(trace.nodes,1), 1]) & ...
    trace.nodes(:,1:3) >= repmat(posM' - sizeM, [size(trace.nodes,1), 1]);
isWithin = isWithin(:,1) & isWithin(:,2) & isWithin(:,3);
if sum(isWithin) > 4
        trace.nodes(~isWithin, :) = 0;
        nodeNrToPlot = find(trace.nodes(:,4));
        edgeNrToPlot = [];
        for i=1:length(nodeNrToPlot)
            edgeNr{i} = find(trace.edges == nodeNrToPlot(i));
            for j=1:length(edgeNr{i})
                if edgeNr{i}(j) > size(trace.edges,1)
                    if sum(trace.edges(edgeNr{i}(j) - size(trace.edges,1),1) == nodeNrToPlot)
                        edgeNrToPlot(end+1) = edgeNr{i}(j) - size(trace.edges,1);
                    end
                else
                    if sum(trace.edges(edgeNr{i}(j),2) == nodeNrToPlot)
                        edgeNrToPlot(end+1) = edgeNr{i}(j);
                    end
                end
            end
        end
        edgeNrToPlot = unique(edgeNrToPlot);
        trace.nodes = trace.nodes(:,1:3) - repmat(posM' - sizeM, [size(trace.nodes, 1) 1]);
        hold on;
        for i=1:length(edgeNrToPlot)
            plot3([trace.nodes(trace.edges(edgeNrToPlot(i),1),2) trace.nodes(trace.edges(edgeNrToPlot(i),2),2)], ...
                [trace.nodes(trace.edges(edgeNrToPlot(i),1),1) trace.nodes(trace.edges(edgeNrToPlot(i),2),1)], ...
                [trace.nodes(trace.edges(edgeNrToPlot(i),1),3) trace.nodes(trace.edges(edgeNrToPlot(i),2),3)], colorString, 'LineWidth', 5);
        end
end
end