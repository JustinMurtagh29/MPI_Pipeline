skeletonFolder = '/gaba/u/kboerg/code/ex145-conversion/';
skeletonFiles = dir(fullfile(skeletonFolder, 'asdf*.nml'));
% define voxel size here as scale in skeletons seems to be wrong (at least in first example I looked at)
voxelSize = [11.24 11.24 28];
result = table;

for i=1:length(skeletonFiles)
    % This file keeps parseNml working for at least few hours (tried twice), check what is happening
    display(skeletonFiles(i).name);
    if ~strcmp(skeletonFiles(i).name, 'asdf82.nml') && ~strcmp(skeletonFiles(i).name, 'asdf83.nml')
        skel = skeleton(fullfile(skeletonFolder, skeletonFiles(i).name));
        skel = skel.deleteEmptyTrees();
        if skel.numTrees == 1
            nodeDegree = skel.calculateNodeDegree();
            assert(all(nodeDegree{1} > 0));
            % Find source (degree 1) and target (degree > 2) nodes
            sourceIdx = find(nodeDegree{1} == 1);
            targetIdx = find(nodeDegree{1} > 2);
            % Find shortest path from all degree one nodes to degree larger 2 nodes
            G = graph(skel.edges{1}(:,1), skel.edges{1}(:,2));
            % Only continue if graph is one connected component
            if all(conncomp(G) == 1)
                d = distances(G, sourceIdx, targetIdx);
                [minDist, closestTargetIdx] = min(d,[],2);
                paths = arrayfun(@(x,y)shortestpath(G,x,y), sourceIdx, targetIdx(closestTargetIdx), 'uni', 0);
                % Determine path length
                pathLength = cellfun(@(x)sum(sqrt(sum(diff(bsxfun(@times, skel.nodes{1}(x,1:3), voxelSize),1,1).^2,2))), paths);
                % Remove all nodes from skeleton that are part of a path which is shorter 5000 nm, note last node (target)
                assert(all(cellfun(@(x)ismember(x(end), targetIdx), paths)));
                pathsWithoutTarget = cellfun(@(x)x(1:end-1), paths, 'uni', 0);
                nodeIdxToRemove = cat(2, pathsWithoutTarget{pathLength < 5000});
                skel = skel.deleteNodes(1, nodeIdxToRemove);
                % Write skeleton to disk and store path length and file length in table
                skel.write(fullfile('/tmpscratch/mberning/gtDendritesSpinesPruned/', skeletonFiles(i).name));
                table{i,2} = skel.physicalPathLength(skel.nodes{1}, skel.edges{1}, voxelSize);
            else
                table{i,3} = 'tree is more than one connected component';
            end
        else
            table{i,3} = 'more than one nonempty tree';
        end
    else
        table{i,3} = 'skeleton not useable with parseNml on gaba?';
    end
    table{i,1} = skeletonFiles(i).name;
end

