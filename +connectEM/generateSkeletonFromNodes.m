function generateSkeletonFromNodes(filename, nodes, treeNames, comments, linearNodes, edges,thr)

    if nargin < 5
        linearNodes = false;
    end
    if ~exist('thr','var')
        thr = 5000; % 5 um edge threshold
    end

    c = 1; % counter for trees in skeleton
    nodeId = 1; % counter for nodes in current skeleton
    nodeOffsetThisSkel = 0;
    % Set colors to be used
    colors = distinguishable_colors(length(nodes), [0 0 0; 1 1 1]);
    colors(:,4) = 0;
    % Write skeleton for check
    skel = initializeSkeleton();
    for tr=1:length(nodes)
        skel{c}.thingID = c;
        skel{c}.name = treeNames{tr};
        skel{c}.color = colors(c,:);
        for no=1:size(nodes{tr},1)
            % Generate comment for skeleton based on agglomeration
            if exist('comments', 'var') && length(comments) >= tr && length(comments{tr}) >= no && ~isempty(comments{tr}{no})
                skel{c}.nodesAsStruct(nodeId-nodeOffsetThisSkel) = generateNodeAsStruct(nodeId, nodes{tr}(no,:), 10, comments{tr}{no});
            else
                skel{c}.nodesAsStruct(nodeId-nodeOffsetThisSkel) = generateNodeAsStruct(nodeId, nodes{tr}(no,:), 10);
            end
            skel{c}.nodes(nodeId-nodeOffsetThisSkel,:) = [nodes{tr}(no,:) 10];
            nodeId = nodeId + 1;
        end
        if exist('edges', 'var')
            skel{c}.edges = edges{tr};
        elseif linearNodes
            skel{c}.edges = cat(1, 1:size(skel{c}.nodes, 1)-1, 2:size(skel{c}.nodes))';
        else
            skel{c}.edges = Graph.getMST(bsxfun(@times, nodes{tr}, [11.24 11.24 28]),thr); 
        end
        c = c + 1;
        nodeOffsetThisSkel = nodeId - 1;
    end

    writeNml(filename, skel, 1);

end

function skel = initializeSkeleton()
    % Set parameters
    skel{1}.parameters.experiment.name='2012-09-28_ex145_07x2_ROI2017';
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


