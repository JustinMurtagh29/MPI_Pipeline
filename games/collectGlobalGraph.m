function collectGlobalGraph(p)

        tic;
        % Collect all edges of local supervoxel graph
        graph.edges = NaN(1e8,2); % 1e8 heuristic for 07x2, turns out actually ~7.5e7
        graph.prob = NaN(1e8,1);
        graph.cubeLI = NaN(1e8,1); % store linear indices as well (second column for correspondences)
        idx = 1;
        for i=1:size(p.local,1)
            for j=1:size(p.local,2)
                for k=1:size(p.local,3)
                    load(p.local(i,j,k).edgeFile);
                    load(p.local(i,j,k).probFile); 
                    % Put everything in one structure
                    nrElements = length(prob);
                    graph.edges(idx:idx+nrElements-1,:) = edges;
                    graph.prob(idx:idx+nrElements-1) = prob;
                    graph.cubeLI(idx:idx+nrElements-1) = repmat(sub2ind(size(p.local),i,j,k), [nrElements 1]);
                    idx = idx+nrElements;
                end
            end
        end

        % Collect all correspondences between local supervoxel graph
        files = dir([p.correspondence.saveFolder, '*global.mat']);
        for i=1:length(files)
            load([p.correspondence.saveFolder files(i).name]);
            nrElements = length(corrGlobal1);
            graph.edges(idx:idx+nrElements-1,:) = [corrGlobal1 corrGlobal2];
            graph.prob(idx:idx+nrElements-1) = ones(length(corrGlobal1),1);
            idx = idx+nrElements;
        end
        % Drop part of arrays preallocated but not assigned
        graph.edges(idx:end,:) = [];
        graph.prob(idx:end) = [];
        graph.cubeLI(idx:end) = [];
        save([p.saveFolder 'graph.mat'], 'graph');
        toc;

end

