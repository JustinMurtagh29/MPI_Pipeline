function generateDendriteQueries(graph, borderMeta, p,segmentMeta)
    %runDirectionality(graph);
    createTasks(borderMeta, p,segmentMeta, graph);
    
end
function runDirectionality(graph)
    options.latentScore = 0.7;
    options.segDirScore = 0.9;
    options.neuriCScore = 0.7;
    options.borderSize = 400;
    options.axonScore = 0;
    options.dendriteScore = 0.5;
    options.recursionSteps = 1;
    options.randomSubset = 10000;
    options.doMerge = false;
    options.startAggloPath = '/gaba/scratch/kboerg/gridAgglo03_514fake.mat';
    gridAgglo_03{514}=load('/gaba/scratch/mberning/aggloGridSearch/search03_00514.mat')
    axonsFinal = gridAgglo_03{514}.dendritesFinal;
    save('/gaba/scratch/kboerg/gridAgglo03_514fake.mat', 'axonsFinal');
    mkdir ('/gaba/scratch/kboerg/dendriteQueries/');
    connectEM.agglomerateDirectionalitySuper2(options, '/gaba/scratch/kboerg/dendriteQueries/', graph)
end
function createTasks(borderMeta, p, segmentMeta, graph)
    direction_col = [];
    startpoint_col = [];
    memory_col = [];
    persistent directionsMB3;
    if isempty(directionsMB3)
        directionsMB3 = load('/gaba/scratch/kboerg/dendriteQueries/1.mat');
    end
    filter = ~cellfun('isempty', directionsMB3.result.latent);
    directionsMB3.result.latent = directionsMB3.result.latent(filter);
    directionsMB3.result.scores = directionsMB3.result.scores(filter);
    directionsMB3.result.pca = directionsMB3.result.pca(filter);
    directionsMB3.result.borderIdx = directionsMB3.result.borderIdx(filter);
    directionsMB3.axons = directionsMB3.axons(filter);
    rng(20170525);
    for aggloidx = randperm(sum(filter));
        if sum(segmentMeta.voxelCount(directionsMB3.axons{aggloidx})) < 200000
            continue;
            %connectEM.skeletonFromAgglo(graph.edges,segmentMeta,directionsMB3.axons(aggloidx),['big',num2str(aggloidx),'_'],'/gaba/scratch/kboerg/dendriteQueries/')
        end
        nodes= useVolume(segmentMeta, directionsMB3, aggloidx);
        todo = directionsMB3.result.latent{aggloidx}(:, 1) > 0.7;
        todo = todo & abs(directionsMB3.result.scores{aggloidx}) > 0.9;
        %  if max(pdist(bsxfun(@times, double(borderMeta.borderCoM(directionsMB3.result.borderIdx{aggloidx}(todo), :)), [11.24, 11.24, 28]))) < 5000
        %      continue;
        %  end
         todo = todo & directionsMB3.result.scores{aggloidx} .* squeeze(sign(directionsMB3.result.pca{aggloidx}(3, 1, :))) > 0.9;
         todoF = find(todo);
         [~, I] = sort(directionsMB3.result.scores{aggloidx}(todoF) .* squeeze(sign(directionsMB3.result.pca{aggloidx}(3, 1, todoF))));
         todoF = todoF(I);
         for endingidx = 1 : length(nodes)
             for idx = 1 : length(todoF)
                 scalize = @(x)bsxfun(@times, x, [11.24,11.24, 28]);
                 if pdist([scalize(nodes(endingidx, :)); scalize(double(borderMeta.borderCoM(directionsMB3.result.borderIdx{aggloidx}(todoF(idx)),:)))]) > 2000
                     continue;
                 end 
                 tic
                 direction = 100 * sign(directionsMB3.result.scores{aggloidx}(todoF(idx))) * directionsMB3.result.pca{aggloidx}(:, 1, todoF(idx))' ./ [11.24, 11.24, 28]; %direction is normalized to 100nm
                 aggloidx
                 startpoint = correctQueryLocationToEndOfSegment(p, directionsMB3.axons{aggloidx}, double(borderMeta.borderCoM(directionsMB3.result.borderIdx{aggloidx}(todoF(idx)),:)), direction, 200, struct('failGracefully', true));
                 bboxSmall = p.bbox + round(bsxfun(@times, 1./[11.24;11.24;28], [3000 -3000; 3000 -3000; 3000 -3000]));
                 if ~(all(bsxfun(@gt, startpoint, bboxSmall(:, 1)'), 2) & all(bsxfun(@lt, startpoint, bboxSmall(:, 2)'), 2))
                     continue;
                 end
                 if double(borderMeta.borderSize(directionsMB3.result.borderIdx{aggloidx}(todoF(idx)))) < 400
                     continue;
                 end
                 startpoint_col(end + 1, :) = startpoint

                 direction_col(end + 1, :) = direction;
                 memory_col(end+1, :) = [aggloidx, idx];
                 toc
                 break
             end
         end
     end
     q.pos={startpoint_col};
     q.dir={direction_col};
     connectEM.generateQueriesFromData(p, segmentMeta, q, '/gaba/scratch/kboerg/testFlightQuery20170524/', struct('writeTasksToFile', true, 'queryBoundingBoxSize', 2000))
end
function nodes = useVolume(segmentMeta, directionsMB3, aggloidx)
    min_coord = connectEM.findQueries(directionsMB3.axons{aggloidx}, segmentMeta);
    load('voxelRepresentation2');
    %addpath('/gaba/u/kboerg/code/skeleton3d-matlab/')

    skel = Skeleton3D(voxelRepresentation2);
    connectEM.querySkeleton
    nodes = [];
    for idx = 1 : length(node2)
        if length(node2(idx).links) == 1
            nodes(end + 1, :) = round([node2(idx).comx * 8, node2(idx).comy*8, node2(idx).comz*8/2.313]+min_coord);
        end
    end
end