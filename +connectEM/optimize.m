% Load parameter struct, graph & segment information needed
load /gaba/u/mberning/results/pipeline/20170217_ROI/allParameterWithSynapses.mat;
graph = load([p.saveFolder 'graph.mat'], 'edges', 'prob');
segMeta = load([p.saveFolder 'segmentMeta.mat'], 'point', 'voxelCount', 'maxSegId');
segMeta.point = segMeta.point';

% Keep only segments larger than 100 voxel
segmentsToKeep = find(segMeta.voxelCount > 100);
idx = all(ismember(graph.edges, segmentsToKeep),2);
% Also keep continuos numbering, should make steps below easier for now
graphCut.edges = connectEM.changem(graph.edges(idx,:), 1:numel(segmentsToKeep), segmentsToKeep);
graphCut.prob = graph.prob(idx);
maxSegId = numel(segmentsToKeep);
clear idx segmentsToKeep;

% Generate better representation
probM = accumarray(double(graphCut.edges), graphCut.prob, ...
    double(repmat(maxSegId, 1, 2)), @max, 0, true);
clear graph*;

threshold = 0.999:-0.001:0.980;
for i=1:length(threshold)
    tic;
    % Inital guess
    x0 = probM > threshold(i);
    [value(i), ~, intraCostPos(i), intraCostNeg(i), ...
        interCostPos(i), interCostNeg(i)] = ...
        connectEM.aggloObjective(probM, x0);
    toc;
end

% Evaluate objective function at starting conditions
% tic; value = connectEM.aggloObjective(probM, x0); toc;
% % How many edges down to X percent still left?
% temp = probM;
% temp(x0) = 0;
% nnz(temp > 0.9)
% There are about 3 million in x0 at 99% (not this does not include CC
% calculated in function) and about 11 million edges left

%% Some first test, switches on edges one at a time
nrIterations = 10000;

xNew = x0;
% Make sure we do not consume too much unnecessary memory
clear x0;
for i=1:nrIterations
    tic;
    % Add random edge with high probability
    if i > 1
        canidates = probM(probM > 0.9);
        canidates(x) = 0;
        [r,c] = find(canidates);
        idx = randi(numel(r), 1);
        xNew(r(idx),c(idx)) = true;
    end
    [value(i), xNew, intraCostPos(i), intraCostNeg(i), ...
        interCostPos(i), interCostNeg(i)] = ...
        connectEM.aggloObjective(probM, xNew);
    % Debug output
    display(['Iteration: ' num2str(i), ', Value: ' num2str(value(i))]);
    display(['IntraPos: ' num2str(intraCostPos(i)), ', IntraNeg: ' num2str(intraCostNeg(i))]);
    display(['InterPos: ' num2str(interCostPos(i)), ', InterNeg: ' num2str(interCostNeg(i))]);
    % Decide whether to reject
    if i > 1
        if value(i-1) > value(i)
            rejected(i) = true;
            display('--> REJECTED');
        else
            rejected(i) = false;
            x = xNew;
        end
    end
    clear xNew;
    toc;
    if mod(i,10) == 0
        display('Saving:');
        tic;
        save(['/gaba/u/mberning/results/aggloObjTests/' num2str(i, '%.5i') '.mat']);
        toc;
    end
end

% %% Try some global optimization approaches (LATER!)
% opts = optimoptions(@fmincon,'Algorithm','interior-point');
% problem = createOptimProblem('fmincon', whatelse);
options = saoptimset('PlotFcns',{@saplotbestx,...
          @saplotbestf,@saplotx,@saplotf});
x = simulannealbnd(fun,x0);

