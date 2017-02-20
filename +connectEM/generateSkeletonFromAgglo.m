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
            % Downsample to 100.000 nodes at most (and use minimal spanning tree)
            %if size(theseCoM, 1) > 10000
            %    idx = randperm(size(theseCoM,1), 10000);
            %    theseCoM = theseCoM(idx,:);
            %    theseEdgesNodes = minimalSpanningTree(theseCoM);
            %    skel{1}.name = [skel{1}.name '_downsampled'];
            %    clear idx;
            %end
            % Write to structure for writeNml
            skel{1}.nodesNumDataAll = zeros(size(theseCoM,1),14);
            skel{1}.nodesNumDataAll(:,1) = 1:size(theseCoM,1); 
            skel{1}.nodesNumDataAll(:,2) = 10*ones(size(theseCoM,1),1);
            skel{1}.nodesNumDataAll(:,3:5) = theseCoM;
            skel{1}.edges = theseEdgesNodes;
            clear theseCoM theseEdgesNodes;
            writeNml([outputFolder treeNames{tr} '.nml'], skel);
            clear skel;
            mappingFile = [outputFolder treeNames{tr} '.txt'];
            script = WK.makeMappingScript(maxSegId, num2cell(cc{tr}));
            fileHandle = fopen(mappingFile, 'w');
            fwrite(fileHandle, script);
            fclose(fileHandle);
        end
    end
end

function skel = initializeSkeleton()
    % Set parameters
    skel{1}.parameters.experiment.name='2012-09-28_ex145_07x2_ROI2016_corrected';
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

function edges = minimalSpanningTree(com)
    if size(com,1) < 2
        edges = [];
    else
        % Minimal spanning tree
        adj = squareform(pdist(bsxfun(@times, com, [11.24 11.24 28])));
        adj(adj > 7200) = 0;
        tree = graphminspantree(sparse(adj), 'Method', 'Kruskal');
        [edges(:,1), edges(:,2)] = find(tree);
    end
end

