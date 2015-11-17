function skel = writeSkeleton(graph, seeds, com);
    c = 1; % counter for trees in skeleton
    nodeId = 1; % counter for nodes in current skeleton
    nodeOffsetThisSkel = 0;
    % Reduce graph only to edges needed here for efficeny
    idx = ismember(graph.edges,cat(1,seeds{:}));
    idx = idx(:,1) & idx(:,2);
    graph.edges = graph.edges(idx,:);
    graph.prob = graph.prob(idx);
    % Write skeleton for check
    skel = initializeSkeleton();
    tic;
    for tr=1:length(seeds)
        if ~isempty(seeds{tr})
            skel{c}.thingID = c;
            skel{c}.name = ['Component ' num2str(tr, '%.2i')];
            skel{c}.color = [rand(1,3) 1];
            theseCoM = com(seeds{tr},:);
            idx = ismember(graph.edges,seeds{tr});
            idx = idx(:,1) & idx(:,2);
            theseEdges = graph.edges(idx,:);
            theseProb = graph.prob(idx);
            for no=1:size(theseCoM,1)
                skel{c}.nodesAsStruct(nodeId-nodeOffsetThisSkel) = generateNodeAsStruct(nodeId, theseCoM(no,:), 10);
                skel{c}.nodes(nodeId-nodeOffsetThisSkel,:) = [theseCoM(no,:) 10];
                nodeId = nodeId + 1;
                theseEdges(theseEdges == seeds{tr}(no)) = no;
            end
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

