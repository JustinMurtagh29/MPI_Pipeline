function addFeatures(parameter, sub)

% get neighbors and connection probabilities....use supervoxelGraph when finished...

load(parameter.local(sub(1),sub(2),sub(3)).segmentFile); % save ids seperately to avoid loading of segments structure!
ids = cell2mat{segments.id};

load(parameter.local(sub(1),sub(2),sub(3)).edgeFile,'edges');
load(parameter.local(sub(1),sub(2),sub(3)).probFile,'prob');

graphData = [edges prob];

for i = 1:length(ids)
    [rows,cols] = ind2sub(size(graphData),find(graphData == ids(i)));
    localEdges = graphData(rows,:);
    if ~isempty(localEdges)
        for j = 1:size(localEdges,1)
            if localEdges(j,1) == ids(i)
                neighbors{i}(j,1) = localEdges(j,2);
            else
                neighbors{i}(j,1) = localEdges(j,1);
            end
            neighbors{i}(j,2) = localEdges(j,3);    
        end
    else
        neighbors{i} = [];
    end
end
    
% filter featuers according to ard and add neighbor features
    % run once with all features to get ard values!
    load(parameter.glia.hyperParameterAllFeatures,'hyp');
    featureFilter = exp(hyp.cov(1:end-1)) < 5; % check cutoff value..
    
    % last analysis Benjamin: including only first max neighbors features gives maximal benefit
    % -> numMax = 1
    numMax = 1;
    
    load(parameter.local(sub(1),sub(2),sub(3)).segmentWeights);
    segmentWeights = segmentWeights(:,featureFilter);
    for i = 1:length(segmentWeights)
        if ~isempty(neighbors{i})
            
            [~,maxIdx(1)] = max(neighbors{i}(:,2));
            neighbors{i}(maxIdx(1),2) = 0;
            [~,maxIdx(2)] = max(neighbors{i}(:,2));
            neighbors{i}(maxIdx(2),2) = 0;
            [~,maxIdx(3)] = max(neighbors{i}(:,2));

            maxId(1) = neighbors{i}(maxIdx(1),1);
            maxId(2) = neighbors{i}(maxIdx(2),1);
            maxId(3) = neighbors{i}(maxIdx(3),1);

            maxSegIdx(1) = find(ids == maxId(1));
            maxSegIdx(2) = find(ids == maxId(2));
            maxSegIdx(3) = find(ids == maxId(3));
             
            if numMax == 1
                weightsNew(i,:) = [segmentWeights(i,:) segmentWeights(maxSegIdx(1),:)];
            elseif numMax == 2
                weightsNew(i,:) = [segmentWeights(i,:) segmentWeights(maxSegIdx(1),:) segmentWeights(maxSegIdx(2),:)];
            else
                weightsNew(i,:) = [segmentWeights(i,:) segmentWeights(maxSegIdx(1),:) segmentWeights(maxSegIdx(2),:) segmentWeights(maxSegIdx(3),:)];
            end
        else
            weightsNew(i,:) = [segmentWeights(i,:) zeros(1,numMax*size(segmentWeights,2))]; 
        end
            
    end
    segmentWeights = weightsNew;
    % introduce new file to save to!!
    save(parameter.local(sub(1),sub(2),sub(3)).segmentWeights, 'segmentWeights');



