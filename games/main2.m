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
            % Cluster (choose too many in order to rather split than merge)
            % Idea is that each process should have approximately 2 border neurites
            % But oversegmented, so it should be overestimated by segIdsB
            % This is time critical step (calculation of eigenvalues)
            tic;
            partition = partitionGraph(graphR, length(segIdsB)./2);
            toc;
            % Write to skeleton file for looking at intermediate result
            skel = writeSkeleton(graph, partition, com);
            writeNml([skeletonSaveDir 'x' num2str(x, '%.2i') 'y' num2str(y, '%.2i') 'z' num2str(z, '%.2i') '.nml'], skel);
            % Find which edges to annotate
            edgesToAnnotate = determineEdgesToAnnotate(graphR, segIdsR, partition);  
            % Write to knowledge database
            v = writeKnowledgeDB(p, graph, edgesToAnnotate);
            % Visualize problems just written in movies

        end
    end
end

