function [output, queryIdx] = detectChiasmataPostSingleNodeLabel(edges, isIntersection, nrExits, nodes, p, nodesV, prob)
temp.edges = edges;
cc = findCCaccordingToGraph(temp, find(isIntersection)); %must be from manuelCode repo
[~, centerOfCC] = cellfun(@(x)min(pdist2(bsxfun(@minus, nodes(x,:), mean(nodes(x,:),1)), [0 0 0])), cc);

% Find out where to query for each CC
queryIdx = cell(length(cc),1);
pos = cell(length(cc),1);
dir = cell(length(cc),1);
if sum(isIntersection) > 100000
    % save('/tmpscratch/kboerg/visX19_0/visX19_1/prep_singlenodelabel.mat','-v7.3');
    % 
    % functionH = @connectEM.detectChiasmataPostSingleNodeLabelSub;
    % inputCell = cellfun(@(x){x}, num2cell(1 : 500), 'uni', 0);
    % cluster = Cluster.getCluster( ...
    %     '-pe openmp 1', ...
    %     '-p 0', ...
    %     '-l h_vmem=24G', ...
    %     '-l s_rt=23:50:00', ...
    %     '-l h_rt=24:00:00');
    % job = Cluster.startJob( functionH, inputCell, ...
    %     'name', 'chiasmata2', ...
    %     'cluster', cluster);
    % Cluster.waitForJob(job);
    for idx = 1 : 500
        if exist(['/tmpscratch/kboerg/visX19_0/visX19_1/temp_singlenodelabel_' num2str(idx) '.mat'], 'file')
            temp = load(['/tmpscratch/kboerg/visX19_0/visX19_1/temp_singlenodelabel_' num2str(idx)]);
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

