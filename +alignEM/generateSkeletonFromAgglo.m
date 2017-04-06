function generateSkeletonFromAgglo(edges, com, cc, treeNames, outputFolder, maxSegId,skelParameters)
    if nargin < 7 || isempty(skelParameters)
        skelParameters = struct();
    end
    % Set colors to be used
    colors = distinguishable_colors(length(cc), [0 0 0; 1 1 1]);
    colors(:,4) = 0;
    for tr=1:length(cc)
        if ~isempty(cc{tr})
            % Generate parameters for skeleton
            skel = initializeSkeleton(skelParameters);
            skel{1}.thingID = 1;
            skel{1}.name = [treeNames{tr} '_' num2str(size(com,1))];
            skel{1}.color = colors(tr,:);
            %little fix for now
            cc{tr} = intersect(cc{tr},edges{tr});
            % Get information to write to skeleton
            theseEdgesNodes = changem(double(edges{tr}), 1:size(com{tr},2), cc{tr});
            % Downsample to 100.000 nodes at most (and use minimal spanning tree)
            %if size(com{tr}, 1) > 10000
            %    idx = randperm(size(com{tr},1), 10000);
            %    com{tr} = com{tr}(idx,:);
            %    theseEdgesNodes = minimalSpanningTree(com{tr});
            %    skel{1}.name = [skel{1}.name '_downsampled'];
            %    clear idx;
            %end
            % Write to structure for writeNml
            skel{1}.nodesNumDataAll = zeros(size(com{tr},2),14);
            skel{1}.nodesNumDataAll(:,1) = 1:size(com{tr},2); 
            skel{1}.nodesNumDataAll(:,2) = 10*ones(size(com{tr},2),1);
            skel{1}.nodesNumDataAll(:,3:5) = com{tr}';
            skel{1}.edges = theseEdgesNodes;
            writeNml(fullfile(outputFolder,[treeNames{tr} '.nml']), skel, 1);
            clear skel;
            mappingFile = fullfile(outputFolder,[treeNames{tr} '.txt']);
            script = WK.makeMappingScript(maxSegId, num2cell(cc{tr}),0);
            fileHandle = fopen(mappingFile, 'w');
            fwrite(fileHandle, script);
            fclose(fileHandle);
        end
    end
end

function skel = initializeSkeleton(parameters)
    % Set parameters
    if nargin < 1 || isempty(parameters)
       parameters = struct(); 
    end
    if ~isfield(parameters,'experiment')
        parameters.experiment.name='2012-09-28_ex145_07x2_ROI2017';
    end
    if ~isfield(parameters,'scale')
        parameters.scale.x = '11.24';
        parameters.scale.y = '11.24';
        parameters.scale.z = '28';
    end
    if ~isfield(parameters,'offset')
        parameters.offset.x = '0';
        parameters.offset.y = '0';
        parameters.offset.z = '0';
    end
    skel{1}. parameters = parameters;
    skel{1}.commentsString = {'<comments></comments>'};
    skel{1}.branchpointsString = {};
    skel{1}.branchpoints = [];
end

% This function is supposed to be present in Matlab Mapping Toolbox, just
% reimplemented here as I was not able to checkout the toolbox
function newmap = changem(map, newcode, oldcode)

    [~, idx] = ismember(map, oldcode);
    newmap = newcode(idx);

end
