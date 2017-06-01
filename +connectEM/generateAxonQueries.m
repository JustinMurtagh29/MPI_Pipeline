function generateAxonQueries(p, graph, segmentMeta, borderMeta, collection)
    direction_col = [];
    startpoint_col = [];
    memory_col = [];
    
    for aggloidx = 1 : size(collection.axons);
        aggloidx
        %there are many small axons fragments
        if max(pdist(bsxfun(@times, double(borderMeta.borderCoM(collection.result.borderIdx{aggloidx}, :)), [11.24, 11.24, 28]))) < 5000
            continue;
        end
        todo = collection.result.latent{aggloidx}(:, 1) > 0.7;
        todo = todo & abs(collection.result.scores{aggloidx}) > 0.9;
        correctlyFlippedScores = collection.result.scores{aggloidx} .* squeeze(sign(collection.result.pca{aggloidx}(3, 1, :)));
        todo = todo &  correctlyFlippedScores > 0.9;
        todoF = find(todo);
        %sort the PCA endings in order how promising they are
        [~, I] = sort(correctlyFlippedScores(todoF));
        todoF = todoF(I);
        for idx = 1 : length(todoF)
            tic
            direction = 100 * sign(collection.result.scores{aggloidx}(todoF(idx))) * collection.result.pca{aggloidx}(:, 1, todoF(idx))' ./ [11.24, 11.24, 28]; %direction is normalized to 100nm
            aggloidx
            startpoint = correctQueryLocationToEndOfSegment(p, collection.axons{aggloidx}, double(borderMeta.borderCoM(collection.result.borderIdx{aggloidx}(todoF(idx)),:)), direction, 200, struct('failGracefully', true));
            % check for border of dataset
            bboxSmall = p.bbox + round(bsxfun(@times, 1./[11.24;11.24;28], [3000 -3000; 3000 -3000; 3000 -3000]));
            if ~(all(bsxfun(@gt, startpoint, bboxSmall(:, 1)'), 2) & all(bsxfun(@lt, startpoint, bboxSmall(:, 2)'), 2))
                continue;
            end
            % border threshold
            startpoint_col(end + 1, :) = startpoint
            direction_col(end + 1, :) = direction;
            memory_col(end+1, :) = [aggloidx, idx];
            toc
            %one query per agglomeration is enough
            break
        end
    end
    [C,ia,ic] = unique(startpoint_col, 'rows');
    q.pos={C};
    q.dir={direction_col(ia, :)};
    mkdir(['/gaba/scratch/kboerg/testFlightQueryAxon20170601/', num2str(idxmain), filesep]);
    connectEM.generateQueriesFromData(p, segmentMeta, q, ['/gaba/scratch/kboerg/testFlightQueryAxon20170601/', num2str(idxmain), filesep], struct('writeTasksToFile', true, 'queryBoundingBoxSize', 2000))

end
function dummy()
end