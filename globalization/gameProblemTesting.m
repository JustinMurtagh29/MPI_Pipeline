function gameProblemTesting(p, coord, lowerT, upperT, nrIterations)
    % Some testing for annotating single axon or dendrite
    % Input: Usual parameter structure, coordinate for start, lower and upper threshold, number of iteration for lookup
    % Load supervoxel in global IDs (edges & probabilities before joining)
    load([p.saveFolder 'graph.mat']);
    % Look up segment ID at inital location
    problems(1).currentNodes = readKnossosRoi(p.seg.root, p.seg.prefix, [coord' coord'], 'uint32', '', 'raw');
    % Empty initialization of structure for iteration
    problems(1).edges = {};
    problems(1).prob = {};
    problems(1).cubeLI = {};
    problems(1).alreadyVisited = [];
    % Iterate
    for i=1:nrIterations
        problems(i+1) = getProblemEdges(problems(i), graph, lowerT, upperT);
    end
    skel = writeProblemsToSkeleton(p,problems);
end

function problemsOut = getProblemEdges(problemsIn, graph, lowerT, upperT) 
    for p=1:length(problemsIn.currentNodes)
        % Choose only edges from current node, above lower threshold and not yet visited
        idx = any(graph.edges == problemsIn.currentNodes(p),2) & graph.prob > lowerT & ~any(ismember(graph.edges, problemsIn.alreadyVisited),2); 
        edgesToCurrentId = graph.edges(idx,:);
        probToCurrentId = graph.prob(idx);
        % Calculate all neighbours
        linIdxSelfInEdges = edgesToCurrentId(:) == problemsIn.currentNodes(p); 
        currentNodes{p} = edgesToCurrentId(~linIdxSelfInEdges);
        edges{p} = edgesToCurrentId;
        prob{p} = probToCurrentId;
        cubeLI{p} = graph.cubeLI(idx);
    end
    problemsOut.currentNodes = vertcat(currentNodes{:});
    problemsOut.edges = edges;
    problemsOut.prob = prob;
    problemsOut.cubeLI = cubeLI;
    problemsOut.alreadyVisited = vertcat(problemsIn.alreadyVisited, problemsIn.currentNodes);
end

function skel = writeProblemsToSkeleton(p, problems);
    % Write skeleton for check
    skel = initializeSkeleton();
    for p=1:length(problems)
        for n=1:length(currentNodes)
        skel{i}.thingID = i;
        skel{i}.name = ['Iteration ' num2str(i, '%.2i')];
        skel{i}.nodes = [];
        skel{i}.edges = [];
        for j=1:length(problems(i).edges)
            for k=1:size(problems(i).edges{j},1)
                load([p.local(problems(i).cubeLI{j}(k)).saveFolder 'CoM.mat']); % Loads CoMAll and ProbAll
                load(p.local(problems(i).cubeLI{j}(k)).edgeFile);
                load(p.local(problems(i).cubeLI{j}(k)).borderFile);
                idx = all(bsxfun(@eq, edges, problems(i).edges{j}(k,:)),2);
                skel{i}.nodes(end+1,:) = CoMAll(idx,1:3);
                skel{i}.nodes(end+1,:) = CoMAll(idx,4:6);
                skel{i}.edges = edges;
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
    nodesAsStruct{1}.id = num2str(id);
    nodesAsStruct{1}.x = num2str(pos(1));
    nodesAsStruct{1}.y = num2str(pos(2));
    nodesAsStruct{1}.z = num2str(pos(3));
    nodesAsStruct{1}.radius = num2str(radius);
    % If comment is passed, add to structure
    if nargin > 3
        nodesAsStruct{1}.comment = comment;
    else
        nodesAsStruct{1}.comment = '';
    end
    % Initalize user behaviour logging to default values 
    nodesAsStruct{1}.inVp = num2str(0);
    nodesAsStruct{1}.inMag = num2str(0);
    nodesAsStruct{1}.time = num2str(0);
end

