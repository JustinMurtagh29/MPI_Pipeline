function debugAgglomerationWithNewQueries(graph, possibleAxons, p, q, outputFolder)
    % Visualize proofread ground truth and all agglomeration(s) overlapping after merging by queries 

    % Use only agglomerates where at least one query will be written
    possibleAxons = possibleAxons(~q.exclude);
    q = structfun(@(x)x(~q.exclude), q, 'uni', 0);
    
    % Choose 100 random current agglomerates where at least one (of the 2) queries is written
    rng default;
    idx = randperm(length(possibleAxons), 100); 
    load([p.saveFolder 'segToPointMap.mat']);
    for i=1:length(idx)
        % agglomerate 
        theseNodes{1} = segToPointMap(possibleAxons{idx(i)},:); 
        theseNames{1} = 'agglomerate';
        theseComments{1} = {};
        for j=1:2
                % Write queries indicating direction (and comment to paste behing link for flight location)
                theseNodes{end+1}(1,:) = round(q.pos{idx(i)}(j,:));
                theseNodes{end}(2,:) = round(q.pos{idx(i)}(j,:) + 10 * q.dir{idx(i)}(j,:));
                theseNames{end+1} = ['query ' num2str(j)];
                [phi, theta, psi] = calculateEulerAngles(-q.dir{idx(i)}(j,:));
                theseComments{end+1}{2} = ['#' num2str(q.pos{idx(i)}(j,1)) ',' num2str(q.pos{idx(i)}(j,2)) ',' num2str(q.pos{idx(i)}(j,3)) ',' ...
                    '1,2,' num2str(phi) ',' num2str(theta) ',' num2str(psi)];
        end
        thisSkel = generateSkeletonFromNodes(p,theseNodes, theseNames, theseComments); 
        expression = ['writeNml(''' outputFolder 'newQueries' num2str(i, '%.3i') '.nml'', thisSkel)'];
        evalc(expression);
        clear this* these*;
    end

end

