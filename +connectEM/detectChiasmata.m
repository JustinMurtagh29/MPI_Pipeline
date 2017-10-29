function output = detectChiasmata(p, nodesV, edges, visualize, outputFolder )
% Detect chiasmata in skeletons based on marching sphere algorithm
% Nodes should be in voxel, scaled here
% Create output folder if it does not exist
if ~isempty(outputFolder) && ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

% Scale to nm
% Make sure edges are unique
[edges, ~, idxE] = unique(edges, 'rows');

% Some parameter for algorithm
p.sphereRadiusOuter = 10000; % in nm
p.sphereRadiusInner = 1000; % in nm
p.voxelSize = [11.24 11.24 28];
nodes=bsxfun(@times,nodesV,p.voxelSize);
% for each node ("marching sphere" approach to merger detection)
isIntersection = false(size(nodes,1),1);
nrExits = zeros(size(nodes,1),1);
if size(nodes, 1) < 1E6
    if ~isempty(edges)
        for i=1:size(nodes,1)
            curNrExits = connectEM.detectChiasmataNodes( ...
                nodes, edges, ones(size(edges,1), 1), p, i);
            
            if curNrExits > 3
                isIntersection(i) = true;
                nrExits(i) = curNrExits;
            end
        end
    end
else
    save([outputFolder 'prep']);
    functionH = @connectEM.detectChiasmataSub;
    inputCell = cellfun(@(x){x}, num2cell(1 : 5000), 'uni', 0);
    cluster = Cluster.getCluster( ...
        '-pe openmp 1', ...
        '-p 0', ...
        '-l h_vmem=24G', ...
        '-l s_rt=23:50:00', ...
        '-l h_rt=24:00:00');
    job = Cluster.startJob( functionH, inputCell, ...
        'name', 'chiasmata1', ...
        'sharedInputs', {outputFolder},  'sharedInputsLocation', 2, ...
        'cluster', cluster);
    Cluster.waitForJob(job);
    for idx = 1 : 5000
        if exist([outputFolder 'temp_' num2str(idx) '.mat'], 'file')
            temp = load([outputFolder 'temp_' num2str(idx)],'nrExits', 'isIntersection');
            nrExits = temp.nrExits + nrExits;
            isIntersection = temp.isIntersection | isIntersection;
        else
            warning(['skipped ' num2str(idx)]);
        end
    end
end

% Find CC of detected intersections according to graph
[output, queryIdx] = connectEM.detectChiasmataPostSingleNodeLabel(edges, isIntersection, nrExits, nodes, p, nodesV, ones(size(edges,1),1),outputFolder);

if visualize
    % Write result to skletons for control (detection of intersections)
    comments = cell(size(nodesV,1), 1);
    comments(isIntersection) = arrayfun(@(x)strcat('intersection, ', num2str(x), ' exits'), nrExits(isIntersection), 'uni', 0);
    writeNml([outputFolder 'skelDetected.nml'], writeSkeletonFromNodesAndEdges({nodesV}, {edges}, {comments}, {'axon'}, {[0 0 1 1]}));
    % Write result to skletons for control (query visualization)
    clear comments;
    n{1} = nodesV(~isIntersection,:);
    % Keep only nodes & edges not part of any intersection
    edgeOffset = cumsum(accumarray(find(~isIntersection), 1))';
    e{1} = edges(~any(ismember(edges, find(isIntersection)),2),:);
    e{1} = edgeOffset(e{1});
    comments{1} = cell(size(n{1},1), 1);
    t{1} = 'original tree without detected intersections';
    col{1} = [0 0 1 1];
    for idx=1:length(queryIdx)
        if any(queryIdx{idx}==-1)
            save([outputFolder 'ohoh'], idx);
            continue;
        end
        n{idx+1} = cat(1, nodesV(output.ccCenterIdx(idx),:), nodesV(output.queryIdx{idx},:));
        e{idx+1} = cat(2, ones(size(n{idx+1},1)-1,1), (2:size(n{idx+1},1))');
        comments{idx+1} = cell(size(n{idx+1},1), 1);
        comments{idx+1}{1} = 'center node of intersection, edges inidcating queries';
        t{idx+1} = ['intersection ' num2str(idx)];
        col{idx+1} = [1 0 0 1];
    end
    writeNml([outputFolder 'skelQueries.nml'], writeSkeletonFromNodesAndEdges(n, e, comments, t, col));
end

if ~isempty(outputFolder)
    save(fullfile(outputFolder, 'result.mat'), 'output');
end

end
