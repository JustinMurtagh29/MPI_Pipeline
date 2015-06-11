function gameProblemTesting(p, coord, lowerT, upperT, nrIterations, outputFile)
    % Some testing for annotating single axon or dendrite
    % Input: Usual parameter structure, coordinate for start, lower and upper threshold, number of iteration for lookup
    % Load supervoxel in global IDs (edges & probabilities before joining)
    load([p.saveFolder 'graph.mat']);
    % Look up segment ID at inital location
    currentNodes{1} = readKnossosRoi(p.seg.root, p.seg.prefix, [coord' coord'], 'uint32', '', 'raw');
    visitedLastIter{1} = currentNodes{1};
    alreadyVisited{1} = [];
    % Iterate
    for i=1:nrIterations
        nodeId = 1;
        for j=1:length(currentNodes{i})
            problemEdges = getProblemEdges(p, currentNodes{i}(j), alreadyVisited{i}, graph, lowerT);
            if isfield(problemEdges, 'neighbours')
                iteration(i).node(nodeId) = problemEdges;
                nodeId = nodeId + 1;
            end
        end
        if i < nrIterations
            allNodes = [iteration(i).node(:).neighbours];
            allNeighbours = [allNodes(:).id];
            allProbabilities = [iteration(i).node(:).probability]; 
            visitedLastIter{i+1} = allNeighbours; 
            currentNodes{i+1} = unique(allNeighbours(cellfun(@(x)max(x), allProbabilities) >  upperT));
            alreadyVisited{i+1} = [alreadyVisited{i} visitedLastIter{i}];
        end
    end
    skel = writeProblemsToSkeleton(p,iteration);
    % Write to nml for visual inspection
    writeNml(outputFile, skel);
    % Write to mat for generating mission.json
    save(strrep(outputFile, '.nml', '.mat'), 'iteration');
end

function problems = getProblemEdges(p, currentNode, alreadyVisited, graph, lowerT) 
    % Choose only edges from current node, above lower threshold and not yet visited
    idx = any(graph.edges == currentNode,2) & graph.prob > lowerT & ~any(ismember(graph.edges, alreadyVisited),2); 
    [edgesToCurrentId, idxU] = unique(graph.edges(idx,:), 'rows');
    probToCurrentId = graph.prob(idx);
    probToCurrentId = probToCurrentId(idxU);
    cubeLI = graph.cubeLI(idx);
    cubeLI = cubeLI(idxU);
    % Calculate ID  and CoM of current starting node
    problems.self.id = currentNode;
    problems.self.CoM = getCoMFromNodeId(p, currentNode, graph) - [1 1 1];
    % Same for all neighbours
    for i=1:size(edgesToCurrentId,1)
        idxInRow = find(edgesToCurrentId(i,:) ~= currentNode);
        problems.neighbours(i).id = edgesToCurrentId(i,idxInRow);
        problems.neighbours(i).CoM = getCoMFromNodeId(p, problems.neighbours(i).id, graph) - [1 1 1];
    end
    % Equivalent edge and probability lists (border between self.id & neighbours.id(i))
    for i=1:size(edgesToCurrentId,1)
        % CubeLI is NaN for correspondences cause 2 cubes are involved
        if isnan(cubeLI(i))
            % Correspondences
            problems.CoM{i} = NaN;
            problems.probability{i} = 1;
        else
            % Normal, GP classified edge
            load(p.local(cubeLI(i)).edgeFile);
            load(p.local(cubeLI(i)).borderFile);
            load(p.local(cubeLI(i)).probFile);
            cubeOffset = p.local(cubeLI(i)).bboxSmall(:,1);
            idxInLocalEdges = all(bsxfun(@eq, edges, edgesToCurrentId(i,:)),2) & prob > lowerT;
            problems.CoM{i} = bsxfun(@plus, round(vertcat(borders(idxInLocalEdges).Centroid)), cubeOffset'- [1 1 1]);
            problems.probability{i} = vertcat(prob(idxInLocalEdges));
        end 
    end
end

function CoM = getCoMFromNodeId(p, currentNode, graph)
    persistent seg cubeLILoaded;
    cubeLI = graph.cubeLI(find(any(graph.edges == currentNode,2) & ~isnan(graph.cubeLI), 1)); 
    cubeOffset = p.local(cubeLI).bboxSmall(:,1);
    % Small conditionals to make loading of segmentation be a little more efficent
    if isempty(cubeLILoaded)
        seg = readKnossosRoi(p.seg.root, p.seg.prefix, p.local(cubeLI).bboxSmall, 'uint32', '', 'raw');
        cubeLILoaded = cubeLI;
    else
        if cubeLI ~= cubeLILoaded
            seg = readKnossosRoi(p.seg.root, p.seg.prefix, p.local(cubeLI).bboxSmall, 'uint32', '', 'raw');
            cubeLILoaded = cubeLI;
        end
    end
    % Calculate CoM from segmentation
    thisSegment = regionprops(seg == currentNode, {'Area' 'Centroid'});
    % Can be multiple connected components in small bounding box, if so place in largest CC
    [~,idxLargest] = max([thisSegment(:).Area]);
    beforeXYcorrection = round(thisSegment(idxLargest).Centroid);
    CoM = beforeXYcorrection([2 1 3]) + cubeOffset';
end

function skel = writeProblemsToSkeleton(p, problems);
    c = 1; % counter for trees in skeleton
    nodeId = 1; % counter for nodes in current skeleton
    nodeOffsetThisSkel = 0;
    % Write skeleton for check
    skel = initializeSkeleton();
    for pr=1:length(problems)
        for no=1:length(problems(pr).node)
            edgeId = 1; % counter for edges in current skeleton
            if length(problems(pr).node(no).neighbours) > 0
                skel{c}.thingID = c;
                skel{c}.name = ['Iteration ' num2str(pr, '%.2i') ', Node ' num2str(problems(pr).node(no).self.id, '%.4i')];
                skel{c}.color = [rand(1,3) 1];
                skel{c}.nodesAsStruct(nodeId-nodeOffsetThisSkel) = generateNodeAsStruct(nodeId, problems(pr).node(no).self.CoM, 10, 'root node');
                skel{c}.nodes(nodeId-nodeOffsetThisSkel,:) = [problems(pr).node(no).self.CoM 10];
                nodeId = nodeId + 1;
                for i=1:length(problems(pr).node(no).neighbours)
                    % If CoM = NaN, then it is a correspodence and does not have a border CoM
                    if isnan(problems(pr).node(no).CoM{i})
                        % Only add target CoM node 
                        targetNodePos = problems(pr).node(no).neighbours(i).CoM;
                        skel{c}.nodesAsStruct(nodeId-nodeOffsetThisSkel) = generateNodeAsStruct(nodeId, targetNodePos, 10, 'correspondence');
                        skel{c}.nodes(nodeId-nodeOffsetThisSkel,:) = [targetNodePos 10];
                        skel{c}.edges(edgeId,:) = [1 nodeId-nodeOffsetThisSkel];
                        nodeId = nodeId + 1;
                        edgeId = edgeId + 1;
                    else
                        nrBorders = size(problems(pr).node(no).CoM{i},1);
                        for j=1:nrBorders
                            % Add border CoM node including probability comment
                            borderNodePos = problems(pr).node(no).CoM{i}(j,:);
                            comment =  ['Probability: ' num2str(problems(pr).node(no).probability{i}(j))];
                            skel{c}.nodesAsStruct(nodeId-nodeOffsetThisSkel) = generateNodeAsStruct(nodeId, borderNodePos, 10, comment);
                            skel{c}.nodes(nodeId-nodeOffsetThisSkel,:) = [borderNodePos 10];
                            skel{c}.edges(edgeId,:) = [1 nodeId-nodeOffsetThisSkel];
                            nodeId = nodeId + 1;
                            edgeId = edgeId + 1;
                            % Add target node CoM every time (cause loops should not be possible in oxalis (?) but placing multiple nodes at same location
                            targetNodePos = problems(pr).node(no).neighbours(i).CoM;
                            skel{c}.nodesAsStruct(nodeId-nodeOffsetThisSkel) = generateNodeAsStruct(nodeId, targetNodePos, 10, '');
                            skel{c}.nodes(nodeId-nodeOffsetThisSkel,:) = [targetNodePos 10];
                            skel{c}.edges(edgeId,:) = [nodeId-nodeOffsetThisSkel-1 nodeId-nodeOffsetThisSkel];
                            nodeId = nodeId + 1;
                            edgeId = edgeId + 1;
                        end
                    end
                end
                % Increase skeleton counter, save NodeID offset for next skeleton
                c = c + 1;
                nodeOffsetThisSkel = nodeId - 1;
            end
        end
    end
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

