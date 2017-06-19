function generateDendriteQueries(idxmain)
    addpath('/gaba/u/kboerg/code/skeleton3d-matlab/')
    addpath('/gaba/u/kboerg/code/skel2graph3d-matlab/')
    addpath('/gaba/u/kboerg/code/manuelCode/games/');
    addpath('/gaba/u/kboerg/code/RESCOPaaS/auxiliary/eulerAngles/')
    addpath('/gaba/u/kboerg/code/pipeline/')
    addpath(genpath('/gaba/u/kboerg/code/auxiliaryMethods/'));
    createTasks(idxmain);
end
function createTasks(idxmain)
    if exist(['/gaba/scratch/kboerg/testFlightQuery20170528/', num2str(idxmain), filesep], 'dir')
        return;
    end
    graph = load('/gaba/u/mberning/results/pipeline/20170217_ROI/graphNew.mat', 'edges', 'prob', 'borderIdx');
    [graph.neighbours, neighboursIdx] = Graph.edges2Neighbors(graph.edges);
    graph.neighProb = cellfun(@(x)graph.prob(x), neighboursIdx, 'uni', 0);
    graph.neighBorderIdx = cellfun(@(x)graph.borderIdx(x), neighboursIdx, 'uni', 0);
    
    load('/gaba/u/mberning/results/pipeline/20170217_ROI/allParameterWithSynapses.mat');
    segmentMeta = load('/gaba/u/mberning/results/pipeline/20170217_ROI/segmentMeta.mat');
    segmentMeta = connectEM.addSegmentClassInformation(p, segmentMeta);
    segmentMeta.point = segmentMeta.point';
    borderMeta = load('/gaba/u/mberning/results/pipeline/20170217_ROI/globalBorder.mat', 'borderSize', 'borderCoM');
    globalSegmentPCA = load('/gaba/u/mberning/results/pipeline/20170217_ROI/globalSegmentPCA.mat', 'covMat');
    scalize = @(x)bsxfun(@times, x, [11.24,11.24, 28]);
    
    direction_col = [];
    startpoint_col = [];
    memory_col = [];
    
    temp= load(['/gaba/scratch/mberning/aggloGridSearch/search03_00514.mat']);
    collection.dendrites = temp.dendritesFinal;
    clear temp;
    temp = load('/gaba/u/kboerg/results/dendrites/direction_03_00514.mat');
    collection.result = temp;
    clear temp;
    %filtering because not all directions might be calculated by first step (if still in testing phase)
    
    filter = ~cellfun('isempty', collection.result.latent);
    collection.result.latent = collection.result.latent(filter);
    collection.result.scores = collection.result.scores(filter);
    collection.result.pca = collection.result.pca(filter);
    collection.result.borderIdx = collection.result.borderIdx(filter);
    collection.dendrites = collection.dendrites(filter);
    
    for aggloidx = idxmain:1000:size(collection.dendrites);
        aggloidx
        %there are many small dendrite fragments, mostly pieces of spines
        if sum(segmentMeta.voxelCount(collection.dendrites{aggloidx})) < 200000
            continue;
        end
        %get skeletonization ends
        nodes= useVolume(segmentMeta, collection, aggloidx);
        
        todo = collection.result.latent{aggloidx}(:, 1) > 0.7;
        todo = todo & abs(collection.result.scores{aggloidx}) > 0.9;
        correctlyFlippedScores = collection.result.scores{aggloidx} .* squeeze(sign(collection.result.pca{aggloidx}(3, 1, :)));
        todo = todo &  correctlyFlippedScores > 0.9;
        todoF = find(todo);
        %sort the PCA endings in order how promising they are
        [~, I] = sort(correctlyFlippedScores(todoF));
        todoF = todoF(I);
        for endingidx = 1 : size(nodes,1)
            for idx = 1 : length(todoF)
                % only allow PCA endings that are close enough to skeletonization ends
                if pdist([scalize(nodes(endingidx, :)); scalize(double(borderMeta.borderCoM(collection.result.borderIdx{aggloidx}(todoF(idx)),:)))]) > 2000
                    continue;
                end 
                tic
                direction = 100 * sign(collection.result.scores{aggloidx}(todoF(idx))) * collection.result.pca{aggloidx}(:, 1, todoF(idx))' ./ [11.24, 11.24, 28]; %direction is normalized to 100nm
                aggloidx
                startpoint = correctQueryLocationToEndOfSegment(p, collection.dendrites{aggloidx}, double(borderMeta.borderCoM(collection.result.borderIdx{aggloidx}(todoF(idx)),:)), direction, 200, struct('failGracefully', true));
                % check for border of dataset
                bboxSmall = p.bbox + round(bsxfun(@times, 1./[11.24;11.24;28], [3000 -3000; 3000 -3000; 3000 -3000]));
                if ~(all(bsxfun(@gt, startpoint, bboxSmall(:, 1)'), 2) & all(bsxfun(@lt, startpoint, bboxSmall(:, 2)'), 2))
                    continue;
                end
                % border threshold
                if double(borderMeta.borderSize(collection.result.borderIdx{aggloidx}(todoF(idx)))) < 200
                    continue;
                end
                
                startpoint_col(end + 1, :) = startpoint
                direction_col(end + 1, :) = direction;
                memory_col(end+1, :) = [aggloidx, idx];
                toc
                %one query per skeletonization ending is enough
                break
            end
        end
    end
    [C,ia,ic] = unique(startpoint_col, 'rows');
    q.pos={C};
    q.dir={direction_col(ia, :)};
    mkdir(['/gaba/scratch/kboerg/testFlightQuery20170528/', num2str(idxmain), filesep]);
    connectEM.generateQueriesFromData(p, segmentMeta, q, ['/gaba/scratch/kboerg/testFlightQuery20170528/', num2str(idxmain), filesep], struct('writeTasksToFile', true, 'queryBoundingBoxSize', 2000))
end
function nodes = useVolume(segmentMeta, collection, aggloidx)
    % get smoothed BW data (and offset)
    [voxelRepresentation2, min_coord] = connectEM.findQueries(collection.dendrites{aggloidx}, segmentMeta, struct('debug', false));

    %create skeleton from BW data
    skel = Skeleton3D(voxelRepresentation2>0);
    node2 = connectEM.querySkeleton(skel)
    nodes = [];
    for idx = 1 : length(node2)
        if length(node2(idx).links) == 1
            nodes(end + 1, :) = round([node2(idx).comx * 8, node2(idx).comy*8, node2(idx).comz*8/2.313]+min_coord); % move again into dataset coords
        end
    end
end