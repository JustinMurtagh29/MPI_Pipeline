function [output, queryIdx] = detectChiasmataPostSingleNodeLabel(edges, isIntersection, nrExits, nodes, p, nodesV, prob,outputFolder)
if sum(isIntersection) > 1
    N = sum(isIntersection);
    if N > 10000
        distances = pdist(nodes(isIntersection,:));
        connections = [];
        for idx = 1 : 100
            connections = [connections, find(distances(idx:100:end)<2000)];
        end
        clear distances
        connections = sort(connections);
        % pdist lists the distances in a special way (see documentation)
        % the next two lines define the borders of the right column of those replies
        idxs_end = cumsum(N-1:-1:1);
        idxs_start = [1, cumsum(N-1:-1:2)+1];
        % sort the connections within these boundaries
        [B, I] = sort([connections,(idxs_end + 0.5)]); % B = A(I);
        I_reverse(I) = 1:length(I);
        % find the separators again in the sorted list
        separators = find(B ~= round(B));
        % make a lookup to find for each entry in which box it is
        bins = repelem(1:(N-1),diff([1,separators]));
        %with that we get the right side for each connection
        rightside = bins(I_reverse(1 : end-1));
        rightside = rightside(1:length(connections));
        %the left side depends where ipwdt is in its box
        leftside_pre = idxs_start(rightside);
        leftside = connections-leftside_pre+rightside+1;
        [cccount,ccpre] =graphconncomp(sparse(leftside, rightside,ones(1,length(connections)),N,N),'Directed',false);
    else
        [cccount,ccpre] =graphconncomp(sparse(squareform(pdist(nodes(isIntersection,:)))<2000));
    end
    lookup = find(isIntersection);
    cc = arrayfun(@(x){lookup(ccpre==x)},1:cccount);
    [~, centerOfCC] = cellfun(@(x)min(pdist2(bsxfun(@minus, nodes(x,:), mean(nodes(x,:),1)), [0 0 0])), cc);
else 
    if sum(isIntersection) == 1
        cc = {find(isIntersection)};
        centerOfCC = 1;
    else
        cc = {};
        centerOfCC = [];
    end
end
% Find out where to query for each CC
queryIdx = cell(length(cc),1);
pos = cell(length(cc),1);
dir = cell(length(cc),1);
if sum(isIntersection) > 100000
    save([outputFolder 'prep_singlenodelabel.mat'],'-v7.3');
    
    functionH = @connectEM.detectChiasmataPostSingleNodeLabelSub;
    inputCell = cellfun(@(x){x}, num2cell(1 : 500), 'uni', 0);
    cluster = Cluster.getCluster( ...
        '-pe openmp 1', ...
        '-p 0', ...
        '-l h_vmem=24G', ...
        '-l s_rt=23:50:00', ...
        '-l h_rt=24:00:00');
    job = Cluster.startJob( functionH, inputCell, ...
        'name', 'chiasmata2', ...
        'sharedInputs', {outputFolder},  'sharedInputsLocation', 2, ...
        'cluster', cluster);
    Cluster.waitForJob(job);
    for idx = 1 : 500
        if exist([outputFolder 'temp_singlenodelabel_' num2str(idx) '.mat'], 'file')
            temp = load([outputFolder 'temp_singlenodelabel_' num2str(idx)]);
            pos(idx:500:length(cc)) = temp.pos(idx:500:length(cc));
            dir(idx:500:length(cc)) = temp.dir(idx:500:length(cc));
            queryIdx(idx:500:length(cc)) = temp.queryIdx(idx:500:length(cc));
        else
            warning(['skipped ' num2str(idx)]);
        end
    end
else
    for i=1:length(cc)
        [pos, dir, queryIdx] = connectEM.detectChiasmataPostSingleNodeLabelSubSub(nodes,edges,prob,p,cc,centerOfCC,pos,dir,queryIdx,i)
    end
end
% Create an output structure
output.nodes = nodesV;
output.edges = edges;
output.prob = [];
output.isIntersection = isIntersection;
output.nrExits = nrExits;
output.ccNodeIdx = cc;
output.ccCenterIdx = cellfun(@(x,y)x(y), cc, num2cell(centerOfCC));
output.queryIdx = queryIdx;
output.position = pos;
output.direction = dir;

