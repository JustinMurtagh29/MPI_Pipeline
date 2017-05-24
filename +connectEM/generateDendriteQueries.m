function generateDendriteQueries(graph)
    runDirectionality(graph);
    % createTasks;
end
function runDirectionality(graph)
    options.latentScore = 0.7;
    options.segDirScore = 0.9;
    options.neuriCScore = 0.7;
    options.borderSize = 30;
    options.axonScore = 0;
    options.recursionSteps = 1;
    options.randomSubset = 100;
    options.startAggloPath = '/gaba/scratch/kboerg/gridAgglo03_514fake.mat';
    gridAgglo_03{514}=load('/gaba/scratch/mberning/aggloGridSearch/search03_00514.mat')
    axonsFinal = gridAgglo_03{514}.dendritesFinal;
    save('/gaba/scratch/kboerg/gridAgglo03_514fake.mat', 'axonsFinal');
    mkdir ('/gaba/scratch/kboerg/dendriteQueries/');
    connectEM.agglomerateDirectionalitySuper2(options, '/gaba/scratch/kboerg/dendriteQueries/', graph)
end
function createTasks
    for aggloidx = agglorand
         todo = directionsMB3.result.latent{aggloidx}(:, 1) > 0.7;
         todo = todo & abs(directionsMB3.result.scores{aggloidx}) > 0.9;
         if max(pdist(bsxfun(@times, double(borderMeta.borderCoM(directionsMB3.result.borderIdx{aggloidx}(todo), :)), [11.24, 11.24, 28]))) < 5000
             continue;
         end

         todo = todo & directionsMB3.result.scores{aggloidx} .* squeeze(sign(directionsMB3.result.pca{aggloidx}(3, 1, :))) > 0.9;
         todoF = find(todo);
         for idx = 1 : length(todoF)
             tic
             direction = 100 * sign(directionsMB3.result.scores{aggloidx}(todoF(idx))) * directionsMB3.result.pca{aggloidx}(:, 1, todoF(idx))' ./ [11.24, 11.24, 28]; %direction is normalized to 100nm
             aggloidx
             startpoint = correctQueryLocationToEndOfSegment(p, directionsMB3.axons{aggloidx}, double(borderMeta.borderCoM(directionsMB3.result.borderIdx{aggloidx}(todoF(idx)),:)), direction, 200, struct('failGracefully', true));
             bboxSmall = p.bbox + round(bsxfun(@times, 1./[11.24;11.24;28], [3000 -3000; 3000 -3000; 3000 -3000]));
             if ~(all(bsxfun(@gt, startpoint, bboxSmall(:, 1)'), 2) & all(bsxfun(@lt, startpoint, bboxSmall(:, 2)'), 2))
                 continue;
             end
             startpoint_col(end + 1, :) = startpoint;

             direction_col(end + 1, :) = direction;
             memory_col(end+1, :) = [aggloidx, idx];
             toc
             break
         end
     end
end
