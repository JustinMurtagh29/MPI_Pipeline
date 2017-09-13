function generateSkeletonFromAggloNew(agglo, treeNames, outputFolder, maxSegId,parameters)
if ~exist(outputFolder,'dir')
    mkdir(outputFolder);
end
if nargin < 7
    parameters = [];
end
if ~exist('maxSegId','var') || isempty(maxSegId)
    maxSegId = max(cat(1,agglo(:).nodes),[],1);
    maxSegId = maxSegId(4);
end
% Set colors to be used
colors = distinguishable_colors(min(length(agglo),100), [0 0 0; 1 1 1]);
colors = repmat(colors, ceil(length(agglo)/100), 1);
colors(:,4) = 0;
for tr=1:length(agglo)
    if ~isempty(agglo(tr).nodes)
        % Generate parameters for skeleton
        theseCoM = agglo(tr).nodes(:,1:3);
        skel = initializeSkeleton(parameters);
        skel{1}.thingID = 1;
        skel{1}.name = [treeNames{tr} '_' num2str(size(theseCoM,1))];
        skel{1}.color = colors(tr,:);
        % Get information to write to skeleton
        
        theseEdgesNodes = agglo(tr).edges;
        % Write to structure for writeNml
        skel{1}.nodesNumDataAll = zeros(size(theseCoM,1),14);
        skel{1}.nodesNumDataAll(:,1) = 1:size(theseCoM,1);
        skel{1}.nodesNumDataAll(:,2) = 10*ones(size(theseCoM,1),1);
        skel{1}.nodesNumDataAll(:,3:5) = theseCoM;
        skel{1}.edges = theseEdgesNodes;
        writeNmlSilent(fullfile(outputFolder,[ treeNames{tr} '.nml']), skel, 1);
        clear skel;
    end
end
mappingFile = fullfile(outputFolder, [treeNames{1} '.txt']);
script = WK.makeMappingScript(maxSegId, connectEM.transformAggloNewOldRepr(agglo), false);
fileHandle = fopen(mappingFile, 'w');
fwrite(fileHandle, script);
fclose(fileHandle);

end

function skel = initializeSkeleton(parameters)
if nargin > 0 && ~isempty(parameters)
    skel{1}.parameters = parameters;
else
    % Set parameters
    skel{1}.parameters.experiment.name='2012-09-28_ex145_07x2_ROI2017';
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
