function generateSkeletonFromAgglo(edges, com, cc, treeNames, outputFolder, maxSegId,parameters)
    if nargin < 7
	parameters = [];
    end
    % Set colors to be used
    colors = distinguishable_colors(min(length(cc),100), [0 0 0; 1 1 1]);
    colors = repmat(colors, ceil(length(cc)/100), 1);
    colors(:,4) = 0;
    for tr=1:length(cc)
        if ~isempty(cc{tr})
            % Generate parameters for skeleton
            skel = initializeSkeleton(parameters);
            skel{1}.thingID = 1;
            skel{1}.name = [treeNames{tr} '_' num2str(size(com,1))];
            skel{1}.color = colors(tr,:);
            % Get information to write to skeleton
            theseCoM = com(cc{tr},:);
            idx = all(ismember(edges,cc{tr}),2);
            theseEdgesSegId = edges(idx,:);
            theseEdgesNodes = changem(double(theseEdgesSegId), 1:size(theseCoM,1), cc{tr});
            clear idx theseEdgesSegId;
           % Write to structure for writeNml
            skel{1}.nodesNumDataAll = zeros(size(theseCoM,1),14);
            skel{1}.nodesNumDataAll(:,1) = 1:size(theseCoM,1);
            skel{1}.nodesNumDataAll(:,2) = 10*ones(size(theseCoM,1),1);
            skel{1}.nodesNumDataAll(:,3:5) = theseCoM;
            skel{1}.edges = theseEdgesNodes;
            clear theseCoM theseEdgesNodes;
            writeNmlSilent(fullfile(outputFolder,[ treeNames{tr} '.nml']), skel, 1);
            clear skel;
        end
    end
    mappingFile = fullfile(outputFolder,[treeNames{1} '.txt']);
    script = WK.makeMappingScript(maxSegId, cc, false);
    fileHandle = fopen(mappingFile, 'w');
    fwrite(fileHandle, script);
    fclose(fileHandle);

end

function skel = initializeSkeleton(parameters)
    if nargin > 0 && ~isempty(parameters)
	skel{1}.parameters = parameters;
    else
    % Set parameters
    skel{1}.parameters.experiment.name='H2_3_v2_U1_SubI';
    skel{1}.parameters.scale.x = '11.24';
    skel{1}.parameters.scale.y = '11.24';
    skel{1}.parameters.scale.z = '28';
    skel{1}.parameters.offset.x = '0';
    skel{1}.parameters.offset.y = '0';
    skel{1}.parameters.offset.z = '0';
    end
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
