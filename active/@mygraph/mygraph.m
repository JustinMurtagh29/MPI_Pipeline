classdef mygraph
    % Implements graph with a symetric sparse matrix as storage for
    % adjaceny (this uses more memory, but seems easier to program)
    properties
        edges;
    end
    methods
        function obj = mygraph()
             obj.edges = sparse([]);
        end
        function obj = insertEdges(obj, edge, weight)
            for i=1:size(edge,1)
                obj.edges(edge(i,1), edge(i,2)) = weight(i);
                obj.edges(edge(i,2), edge(i,1)) = weight(i);
            end
        end
        function obj = removeEdges(obj, edge)
            for i=1:size(edge,1)
                obj.edges(edge(i,1), edge(i,2)) = 0;
                obj.edges(edge(i,2), edge(i,1)) = 0;
            end
        end
        function [list, weights] = findNeighbours(obj, nodeNr)
            list = find(obj.edges(nodeNr,:));
            weights = full(obj.edges(nodeNr,list));
            [weights, permutation] = sort(weights);
            list = list(permutation);
        end
        function nrNodes = size(obj)
            nrNodes = size(obj.edges,1);
        end
        function obj = thresholdConnected(obj, thres)
            remain = obj.edges > thres & obj.edges ~= 0;
            newEdges = sparse(size(obj.edges,1), size(obj.edges,2));
            newEdges(remain) = obj.edges(remain);
            singles = find(sum(remain,2) == 0);
            for i=1:length(singles)
                [idx, weights] = findNeighbours(obj, singles(i));
                newEdges(singles(i), idx(end)) = weights(end);
                newEdges(idx(end), singles(i)) = weights(end);
            end
            obj.edges = newEdges;
        end
        function plotGraph(obj, xyz)
            [i,j] = find(obj.edges);
            [~, p] = sort(max(i,j));
            i = i(p);
            j = j(p);
            minimum = min(obj.edges(obj.edges ~= 0));
            range = max(obj.edges(obj.edges ~= 0))-minimum;
            for k=1:length(i)
                X = [ xyz(i(k),1) xyz(j(k),1)]';
                Y = [ xyz(i(k),2) xyz(j(k),2)]';
                Z = [ xyz(i(k),3) xyz(j(k),3)]';
                temp = 0.5+0.5*(obj.edges(i(k),j(k))-minimum)/range;
                C = [1-temp temp 0];
                plot3(X,Y,Z,'LineStyle','-','Color',C,'Marker','*');
                hold on;
            end
        end
    end
end
