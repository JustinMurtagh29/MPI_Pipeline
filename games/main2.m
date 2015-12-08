% Load needed data files
load('/gaba/u/mberning/results/pipeline/20151111T183414/allParameter.mat');
load([p.saveFolder 'graph.mat']);
load([p.saveFolder 'CoM.mat']);
load([p.saveFolder 'agglomeration' filesep 'nucleiVesselBorder.mat']);

skeletonSaveDir = [p.saveFolder 'kdb' filesep 'skeletons' filesep];
if ~exist(skeletonSaveDir, 'dir')
    mkdir(skeletonSaveDir);
end

for x=6%1:size(p.local,1)
    for y=9%1:size(p.local,2)
        for z=9%1:size(p.local,3)
            seedCube = p.local(x,y,z).bboxSmall;
            % Restrict to region & remove somata and vessel indices
            excludeIds = cat(1, agglo.nucleiGrown{:}, agglo.vessel{:});
            [graphR, comR, segIdsR, segIdsB] = restrictSGtoRegion(p, graph, com, seedCube, excludeIds); 
            % Remove all edges below some threshold (will be subdivded by partioning, so lower is ok)
            [graphRT, unusedIds] = restrictSGtoAboveProbability(graphR, .95);
            % Determine number of connected components in thresholded graph
            cc = Graph.findConnectedComponents(graphRT.edges, false, true);
            % Cluster (choose too many in order to rather split than merge)
            % Idea is that each process should have approximately 2 border neurites
            % So it should be overestimated (factor 2) by segIdsB
            % This is time critical step (calculation of eigenvalues & kmeans)
            %if length(cc) > floor(0.6*length(segIdsB))
            %    error('Partioning will not work if #CC > #wantedPartitions'); 
            %end
            %tic;
            %partition = partitionGraph(graphRT, floor(0.6*length(segIdsB)));
            %toc;
            % Seperate partition by CC on SG (this should not happen in different clustering method?)
            %partitionCorrected = splitPartitionByCC(graphRT, partition);
            % Write to skeleton file for looking at intermediate result
            skel = writeSkeleton(graph, cc, com);
            writeNml([skeletonSaveDir 'x' num2str(x, '%.2i') 'y' num2str(y, '%.2i') 'z' num2str(z, '%.2i') '_95perCC.nml'], skel);
            % Find which edges to annotate
            [edgesToAnnotate, probToAnnoate, edgesToForget, probToForget] = determineEdgesToAnnotate(graphR, unusedIds);
            cc1 = Graph.findConnectedComponents(edgesToAnnotate, false, true);
            skel = writeSkeletonEdges(graph, com, cc1, edgesToAnnotate, probToAnnoate, unusedIds);
            writeNml([skeletonSaveDir 'x' num2str(x, '%.2i') 'y' num2str(y, '%.2i') 'z' num2str(z, '%.2i') '_edgesToAnnotate.nml'], skel);
            cc2 = Graph.findConnectedComponents(edgesToForget, false, true);
            skel = writeSkeletonEdges(graph, com, cc2, edgesToForget, probToForget, unusedIds);
            writeNml([skeletonSaveDir 'x' num2str(x, '%.2i') 'y' num2str(y, '%.2i') 'z' num2str(z, '%.2i') '_edgesToForget.nml'], skel);
            % Write to knowledge database
            v = writeKnowledgeDB(p, graph, edgesToAnnotate);
            % Visualize problems just written in movies

        end
    end
end

