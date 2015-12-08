function skel = writeSkeletonEdges(graph, com, cc, probabilities);
      
    c = 1; % counter for trees in skeleton
    nodeId = 1; % counter for nodes in current skeleton
    nodeOffsetThisSkel = 0;
    % Write skeleton for check
    skel = initializeSkeleton();
    tic;
    for tr=1:length(cc)
        if ~isempty(cc{tr})
            skel{c}.thingID = c;
            skel{c}.name = ['Component ' num2str(tr, '%.4i')];
            skel{c}.color = [rand(1,3) 1];
            theseCoM = com(cc{tr},:);
            idx = ismember(graph.edges,cc{tr});
            idx = idx(:,1) & idx(:,2);
            theseEdges = graph.edges(idx,:);
            theseProb = graph.prob(idx);
            for no=1:size(theseCoM,1)
                if no > 1
                    cS = ['Step: ' num2str(no-1, '%.3i') ' Probability: ' num2str(probabilities{tr}(no), '%.2f')];
                    skel{c}.nodesAsStruct(nodeId-nodeOffsetThisSkel) = generateNodeAsStruct(nodeId, theseCoM(no,:), 10, cS);
                else
                    skel{c}.nodesAsStruct(nodeId-nodeOffsetThisSkel) = generateNodeAsStruct(nodeId, theseCoM(no,:), 10);
                end
                skel{c}.nodes(nodeId-nodeOffsetThisSkel,:) = [theseCoM(no,:) 10];
                nodeId = nodeId + 1;
                theseEdges(theseEdges == cc{tr}(no)) = no;
            end
            skel{c}.edges = theseEdges;
            c = c + 1;
            nodeOffsetThisSkel = nodeId - 1;
            Util.progressBar(tr, length(cc));
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

