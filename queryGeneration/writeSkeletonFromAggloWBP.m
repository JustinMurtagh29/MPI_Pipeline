function skel = writeSkeletonFromAggloWBP(graph, com, cc, excluded, seedIds, querryIds, treeName)

    c = 1; % counter for trees in skeleton
    nodeId = 1; % counter for nodes in current skeleton
    nodeOffsetThisSkel = 0;
    % Set colors to be used
    colors = distinguishable_colors(length(cc), [0 0 0; 1 1 1]);
    colors(:,4) = 0;
    % Write skeleton for check
    skel = initializeSkeleton();
    for tr=1:length(cc)
        if ~isempty(cc{tr})
            skel{c}.thingID = c;
            skel{c}.name = treeName;
            skel{c}.color = colors(c,:);
            theseCoM = com(cc{tr},:);
            idx = ismember(graph.edges,cc{tr});
            idx = idx(:,1) & idx(:,2);
            theseEdges = graph.edges(idx,:);
            for no=1:size(theseCoM,1)
                % Save nodes from agglomeration
                idx = find(cc{tr}(no) == querryIds);
                if idx
                    comment = ['querry pos ' num2str(idx)];
                else
                    comment = '';
                end
                skel{c}.nodesAsStruct(nodeId-nodeOffsetThisSkel) = generateNodeAsStruct(nodeId, theseCoM(no,:), 10, comment);
                skel{c}.nodes(nodeId-nodeOffsetThisSkel,:) = [theseCoM(no,:) 10];
                nodeId = nodeId + 1;
                theseEdges(theseEdges == cc{tr}(no)) = no;
            end
            [uniqueIds, idx] = unique(cat(2,(excluded(:).ids)));
            reasons = cat(2,excluded(:).reasons);
            uniqueReasons = reasons(idx);
            for i=1:length(uniqueIds)
                % Also save nodes that have been excluded for now for better debugging, mark with comment
                skel{c}.nodesAsStruct(nodeId-nodeOffsetThisSkel) = generateNodeAsStruct(nodeId, com(uniqueIds(i),:), 10, uniqueReasons{i});
                skel{c}.nodes(nodeId-nodeOffsetThisSkel,:) = [com(uniqueIds(i),:) 10];
                nodeId = nodeId + 1;
                [~, minIdx] = min(pdist2(bsxfun(@times, com(uniqueIds(i),:), [11.24 11.24 28]), bsxfun(@times, theseCoM, [11.24 11.24 28])));
                % Generate comment for skeleton based on agglomeration
                theseEdges(end+1,:) = [minIdx size(theseCoM,1)+i];
            end
            skel{c}.edges = theseEdges; 
            skel{1}.branchpoints = find(ismember(cc{tr}, seedIds)) + nodeOffsetThisSkel;
            c = c + 1;
            nodeOffsetThisSkel = nodeId - 1;
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

