function generateSkeletonFromAgglo(edges, com, cc, treeNames, outputFolder, maxSegId)

    % Set colors to be used
    colors = distinguishable_colors(length(cc), [0 0 0; 1 1 1]);
    colors(:,4) = 0;
    for tr=1:length(cc)
        if ~isempty(cc{tr})
            % Generate parameters for skeleton
            skel = initializeSkeleton();
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
            writeNmlSilent([outputFolder treeNames{tr} '.nml'], skel, 1);
            clear skel;
        end
    end
    mappingFile = [outputFolder treeNames{1} '.txt'];
    script = WK.makeMappingScript(maxSegId, cc, false);
    fileHandle = fopen(mappingFile, 'w');
    fwrite(fileHandle, script);
    fclose(fileHandle);

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

% This function is supposed to be present in Matlab Mapping Toolbox, just
% reimplemented here as I was not able to checkout the toolbox
function newmap = changem(map, newcode, oldcode)

    [~, idx] = ismember(map, oldcode);
    newmap = newcode(idx);

end
