function generateDendriteQueries2(idxmain, graph,segmentMeta,borderMeta,globalSegmentPCA)
    addpath('/gaba/u/kboerg/code/skeleton3d-matlab/')
    addpath('/gaba/u/kboerg/code/skel2graph3d-matlab/')
    addpath('/gaba/u/kboerg/code/manuelCode/games/');
    addpath('/gaba/u/kboerg/code/RESCOPaaS/auxiliary/eulerAngles/')
    addpath('/gaba/u/kboerg/code/pipeline/')
    addpath(genpath('/gaba/u/kboerg/code/auxiliaryMethods/'));
    %runDirectionality(graph);
    %createTasks(borderMeta, p,segmentMeta, graph, agglos);
    createTasks(idxmain,graph,segmentMeta,borderMeta,globalSegmentPCA);
    %y = connectEM.calculateDirectionalityOfAgglomerates(agglo{end}, graph, segmentMeta, borderMeta, globalSegmentPCA, struct('bboxDist',2000, 'voxelSize', [11.24, 11.24, 28]));
    %save('tempy', 'y', 'v);
end
function runDirectionality(graph)
    options.latentScore = 0.7;
    options.segDirScore = 0.9;
    options.neuriCScore = 0.7;
    options.borderSize = 400;
    options.axonScore = 0;
    options.dendriteScore = 0.5;
    options.recursionSteps = 1;
    options.randomSubset = 100;
    options.doMerge = false;
    options.startAggloPath = '/gaba/scratch/mberning/aggloGridSearch/search03_00514.mat';
    %gridAgglo_03{514}=load('/gaba/scratch/mberning/aggloGridSearch/search03_00514.mat')
    %axonsFinal = gridAgglo_03{514}.dendritesFinal;
    %save('/gaba/scratch/kboerg/gridAgglo03_514fake.mat', 'axonsFinal');
    mkdir ('/gaba/scratch/kboerg/dendriteQueries/');
    connectEM.agglomerateDirectionalitySuper2(options, '/gaba/scratch/kboerg/dendriteQueries/', graph)
end
function createTasks(idxmain,graph,segmentMeta,borderMeta,globalSegmentPCA)
    % if exist(['/gaba/scratch/kboerg/testFlightQuery20170528/', num2str(idxmain), filesep], 'dir')
    %     return;
    % end
    % graph = load('/gaba/u/mberning/results/pipeline/20170217_ROI/graphNew.mat', 'edges', 'prob', 'borderIdx');
    % [graph.neighbours, neighboursIdx] = Graph.edges2Neighbors(graph.edges);
    % graph.neighProb = cellfun(@(x)graph.prob(x), neighboursIdx, 'uni', 0);
    % graph.neighBorderIdx = cellfun(@(x)graph.borderIdx(x), neighboursIdx, 'uni', 0);
    % 
    load('/gaba/u/mberning/results/pipeline/20170217_ROI/allParameterWithSynapses.mat');
    % segmentMeta = load('/gaba/u/mberning/results/pipeline/20170217_ROI/segmentMeta.mat');
    % segmentMeta = connectEM.addSegmentClassInformation(p, segmentMeta);
    % segmentMeta.point = segmentMeta.point';
    % borderMeta = load('/gaba/u/mberning/results/pipeline/20170217_ROI/globalBorder.mat', 'borderSize', 'borderCoM');
    % globalSegmentPCA = load('/gaba/u/mberning/results/pipeline/20170217_ROI/globalSegmentPCA.mat', 'covMat');
    scalize = @(x)bsxfun(@times, x, [11.24,11.24, 28]);
    
    direction_col = [];
    startpoint_col = [];
    memory_col = [];
    %persistent directionsMB3;
    %if isempty(directionsMB3)
    %    directionsMB3 = load('/gaba/scratch/kboerg/dendriteQueries/1.mat');
    %end
    zz= load(['/gaba/scratch/mberning/aggloGridSearch/search03_00514.mat']);
    directionsMB3.axons = zz.dendritesFinal;
    clear zz;
    y=connectEM.calculateDirectionalityOfAgglomerates(directionsMB3.axons(5545),graph,segmentMeta,borderMeta,globalSegmentPCA,struct('bboxDist',2000,'voxelSize',[11.24,11.24,28]));
    directionsMB3.result = y;
    clear y;
    directionsMB3.axons(1)= directionsMB3.axons(5545);
    %filtering because not all directions might be calculated by first step (if still in testing phase)
    
    % filter = ~cellfun('isempty', directionsMB3.result.latent);
    % directionsMB3.result.latent = directionsMB3.result.latent(filter);
    % directionsMB3.result.scores = directionsMB3.result.scores(filter);
    % directionsMB3.result.pca = directionsMB3.result.pca(filter);
    % directionsMB3.result.borderIdx = directionsMB3.result.borderIdx(filter);
    % directionsMB3.axons = directionsMB3.axons(filter);
    %shuffling for testing, so that first piece is more representative
    %rng(20170525);
    for aggloidx = idxmain:1000:100000%randperm(sum(filter));
        aggloidx
        %there are many small dendrite fragments, mostly pieces of spines
        if sum(segmentMeta.voxelCount(directionsMB3.axons{aggloidx})) < 200000
            continue;
            %connectEM.skeletonFromAgglo(graph.edges,segmentMeta,agglo{end}.dendritesFinal(5858),['weirdness',num2str(aggloidx),'_'],'/gaba/scratch/kboerg/dendriteQueries/')
        end
        %get skeletonization ends
        nodes= useVolume(segmentMeta, directionsMB3, aggloidx);
        
        todo = directionsMB3.result.latent{aggloidx}(:, 1) > 0.7;
        todo = todo & abs(directionsMB3.result.scores{aggloidx}) > 0.9;
        correctlyFlippedScores = directionsMB3.result.scores{aggloidx} .* squeeze(sign(directionsMB3.result.pca{aggloidx}(3, 1, :)));
        todo = todo &  correctlyFlippedScores > 0.9;
        todoF = find(todo);
        %sort the PCA endings in order how promising they are
        [~, I] = sort(correctlyFlippedScores(todoF));
        todoF = todoF(I);
        for endingidx = 1 : size(nodes,1)
            for idx = 1 : length(todoF)
                % only allow PCA endings that are close enough to skeletonization ends
                if pdist([scalize(nodes(endingidx, :)); scalize(double(borderMeta.borderCoM(directionsMB3.result.borderIdx{aggloidx}(todoF(idx)),:)))]) > 2000
                    continue;
                end 
                tic
                direction = 100 * sign(directionsMB3.result.scores{aggloidx}(todoF(idx))) * directionsMB3.result.pca{aggloidx}(:, 1, todoF(idx))' ./ [11.24, 11.24, 28]; %direction is normalized to 100nm
                aggloidx
                startpoint = correctQueryLocationToEndOfSegment(p, directionsMB3.axons{aggloidx}, double(borderMeta.borderCoM(directionsMB3.result.borderIdx{aggloidx}(todoF(idx)),:)), direction, 200, struct('failGracefully', true));
                % check for border of dataset
                bboxSmall = p.bbox + round(bsxfun(@times, 1./[11.24;11.24;28], [3000 -3000; 3000 -3000; 3000 -3000]));
                if ~(all(bsxfun(@gt, startpoint, bboxSmall(:, 1)'), 2) & all(bsxfun(@lt, startpoint, bboxSmall(:, 2)'), 2))
                    continue;
                end
                % border threshold
                if double(borderMeta.borderSize(directionsMB3.result.borderIdx{aggloidx}(todoF(idx)))) < 200
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
function nodes = useVolume(segmentMeta, directionsMB3, aggloidx)
    % get smoothed BW data (and offset)
    [voxelRepresentation3, min_coord] = connectEM.findQueries(directionsMB3.axons{aggloidx}, segmentMeta, struct('debug', false));

    %create skeleton from BW data
    skel = Skeleton3D(voxelRepresentation3>0);
    node2 = connectEM.querySkeleton(skel)
    nodes = [];
    for idx = 1 : length(node2)
        if length(node2(idx).links) == 1
            nodes(end + 1, :) = round([node2(idx).comx * 8, node2(idx).comy*8, node2(idx).comz*8/2.313]+min_coord); % move again into dataset coords
        end
    end
end