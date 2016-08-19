function skel = writeSkeletonFromNodes2(p,nodes, treeNames, comments)

    if nargin < 4
        treePrefix = 'Component ';
    end

    c = 1; % counter for trees in skeleton
    nodeId = 1; % counter for nodes in current skeleton
    nodeOffsetThisSkel = 0;
    % Set colors to be used
    colors = distinguishable_colors(length(nodes), [0 0 0; 1 1 1]);
    colors(:,4) = 0;
    % Write skeleton for check
    skel = initializeSkeleton(p);
    for tr=1:length(nodes)
        skel{c}.thingID = c;
        skel{c}.name = treeNames{tr};
        skel{c}.color = colors(c,:);
        for no=1:size(nodes{tr},1)
            % Generate comment for skeleton based on agglomeration
            if length(comments) >= tr && length(comments{tr}) >= no && ~isempty(comments{tr}{no})
                skel{c}.nodesAsStruct(nodeId-nodeOffsetThisSkel) = generateNodeAsStruct(nodeId, nodes{tr}(no,:), 10, comments{tr}{no});
            else
                skel{c}.nodesAsStruct(nodeId-nodeOffsetThisSkel) = generateNodeAsStruct(nodeId, nodes{tr}(no,:), 10);
            end
            skel{c}.nodes(nodeId-nodeOffsetThisSkel,:) = [nodes{tr}(no,:) 10];
            nodeId = nodeId + 1;
        end
        skel{c}.edges = minimalSpanningTree(p,nodes{tr}); 
        c = c + 1;
        nodeOffsetThisSkel = nodeId - 1;
    end
end

function skel = initializeSkeleton(p)
    % Set parameters
    skel{1}.parameters.experiment.name=p.experimentName;
    skel{1}.parameters.scale.x = num2str(p.raw.voxelSize(1));
    skel{1}.parameters.scale.y = num2str(p.raw.voxelSize(2));
    skel{1}.parameters.scale.z = num2str(p.raw.voxelSize(3));
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

function edges = minimalSpanningTree(p,com)
    if size(com,1) < 2
        edges = [];
    else
        % Minimal spanning tree
        adj = squareform(pdist(bsxfun(@times, com, [p.raw.voxelSize(1) p.raw.voxelSize(2) p.raw.voxelSize(3)])));
        tree = graphminspantree(sparse(adj), 'Method', 'Kruskal');
        [edges(:,1), edges(:,2)] = find(tree);
    end
end
