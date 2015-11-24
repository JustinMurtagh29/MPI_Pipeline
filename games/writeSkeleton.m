function skel = writeSkeleton(graph, seeds, com, writeMaxProbableContFlag, startSeeds);
    
    if nargin < 4
        writeMaxProbableContFlag = false;
    end    
    
    c = 1; % counter for trees in skeleton
    nodeId = 1; % counter for nodes in current skeleton
    nodeOffsetThisSkel = 0;
    % Reduce graph only to edges needed here for efficeny
    idx = ismember(graph.edges,cat(1,seeds{:}));
    idx = idx(:,1) & idx(:,2);
    edges = graph.edges(idx,:);
    prob = graph.prob(idx);
    % Write skeleton for check
    skel = initializeSkeleton();
    tic;
    for tr=1:length(seeds)
        if ~isempty(seeds{tr})
            skel{c}.thingID = c;
            skel{c}.name = ['Component ' num2str(tr, '%.2i')];
            skel{c}.color = [rand(1,3) 1];
            theseCoM = com(seeds{tr},:);
            idx = ismember(edges,seeds{tr});
            idx = idx(:,1) & idx(:,2);
            theseEdges = edges(idx,:);
            theseProb = prob(idx);
            for no=1:size(theseCoM,1)
                skel{c}.nodesAsStruct(nodeId-nodeOffsetThisSkel) = generateNodeAsStruct(nodeId, theseCoM(no,:), 10);
                skel{c}.nodes(nodeId-nodeOffsetThisSkel,:) = [theseCoM(no,:) 10];
                nodeId = nodeId + 1;
                theseEdges(theseEdges == seeds{tr}(no)) = no;
            end
            if writeMaxProbableContFlag
                % Find endSegment (if more than 1 segment use PC)
                if length(seeds{tr}) > 1
                    % If onlt two segments PCA will not work
                    if length(seeds{tr}) == 2
                        maxIdx = 1;
                        minIdx = 2;
                    else
                        [~,score] = princomp(theseCoM);
                        [~, maxIdx] = max(score(:,1));
                        [~, minIdx] = min(score(:,1));
                    end
                    if ismember(seeds{tr}(maxIdx), [startSeeds{:}]);
                        endSegment = seeds{tr}(minIdx);
                        edgeIdxInNml = minIdx;
                    elseif ismember(seeds{tr}(minIdx), [startSeeds{:}])
                        endSegment = seeds{tr}(maxIdx);
                        edgeIdxInNml = maxIdx;
                    % These two next elseif only work for current plane
                    elseif com(seeds{tr}(maxIdx),2) > com(seeds{tr}(minIdx),2)
                        endSegment = seeds{tr}(minIdx);
                        edgeIdxInNml = minIdx;
                    elseif com(seeds{tr}(minIdx),2) > com(seeds{tr}(maxIdx),2)
                        endSegment = seeds{tr}(maxIdx);
                        edgeIdxInNml = maxIdx;
                    % Else choose random for now
                    else
                        rN = randi(2, 1);
                        if rN == 1
                            endSegment = seeds{tr}(minIdx);
                            edgeIdxInNml = minIdx;                           
                        else
                            endSegment = seeds{tr}(minIdx);
                            edgeIdxInNml = minIdx;
                        end
                    end
                else
                    endSegment = seeds{tr};
                    edgeIdxInNml = 1;
                end
                % Find maximal probability connection to end segment
                idxToEnd = any(ismember(graph.edges, endSegment),2);
                idxInAgglo = all(ismember(graph.edges, cat(1, seeds{:})),2);
                idx = find(idxToEnd & ~idxInAgglo);
                if isempty(idx)
                    skel{c}.nodesAsStruct{nodeId-1-nodeOffsetThisSkel}.comment = 'no (not already collected) neighbours found';
                else
                    [maxProb,maxProbIdx] = max(graph.prob(idx));
                    querrySegment = graph.edges(idx(maxProbIdx),:);
                    querrySegment = querrySegment(querrySegment ~= endSegment);
                    % Update edge and probability lists (with additional edge to annotate)
                    theseEdges(end+1,:) = [edgeIdxInNml; max(theseEdges(:))+1];
                    theseProb(end+1) = maxProb;
                    % Write node to skeleton
                    comment = ['Destination of querry: edge p=' num2str(maxProb, '%3.2f')];
                    skel{c}.nodesAsStruct(nodeId-nodeOffsetThisSkel) = generateNodeAsStruct(nodeId, com(querrySegment,:), 10, comment);
                    skel{c}.nodes(nodeId-nodeOffsetThisSkel,:) = [com(querrySegment,:) 10];
                    nodeId = nodeId + 1;
                end
            end
            % Generate minimal spanning tree if necessary
            if size(theseEdges,1) < 3
                skel{c}.edges = theseEdges;
            else
                theseEdges = minimalSpanningTree(theseEdges, 2-theseProb); % to keep correspondence cost positive
                skel{c}.edges = theseEdges;
            end
            c = c + 1;
            nodeOffsetThisSkel = nodeId - 1;
            Util.progressBar(tr, length(seeds));
        end
    end
end

function edgesNew = minimalSpanningTree(edges, prob)
    maxID = max(edges(:));
    adj = sparse(edges(:,1), edges(:,2), prob, maxID, maxID);
    adj = adj + adj';
    tree = graphminspantree(adj);
    [edgesNew(:,1), edgesNew(:,2)] = find(tree);
end

function skel = initializeSkeleton()
    % Set parameters
    skel{1}.parameters.experiment.name='2012-09-28_ex145_07x2_segNew';
    skel{1}.parameters.scale.x = '11.24';
    skel{1}.parameters.scale.y = '11.24';
    skel{1}.parameters.scale.z = '28';
    skel{1}.parameters.offset.x = '0';
    skel{1}.parameters.offset.y = '0';
    skel{1}.parameters.offset.z = '0';
    skel{1}.commentsString = {'<comments></comments>'};
    skel{1}.branchpointsString = {};
    skel{1}.branchpoints = [];
end

function nodeAsStruct = generateNodeAsStruct(id,pos,radius,comment)
    % Why would one store everything in strings, well let's do it anyway :)
    % Important stuff, id, position, radius
    nodeAsStruct{1}.id = num2str(id);
    nodeAsStruct{1}.x = num2str(pos(1));
    nodeAsStruct{1}.y = num2str(pos(2));
    nodeAsStruct{1}.z = num2str(pos(3));
    nodeAsStruct{1}.radius = num2str(radius);
    % If comment is passed, add to structure
    if nargin > 3
        nodeAsStruct{1}.comment = comment;
    else
        nodeAsStruct{1}.comment = '';
    end
    % Initalize user behaviour logging to default values 
    nodeAsStruct{1}.inVp = num2str(0);
    nodeAsStruct{1}.inMag = num2str(0);
    nodeAsStruct{1}.time = num2str(0);
end

