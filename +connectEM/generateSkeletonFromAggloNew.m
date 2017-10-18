function allSkel = generateSkeletonFromAggloNew(agglo, treeNames, outputFolder, maxSegId,parameters,filename)
% if filename is given, all skeletons are written into one file

if ~exist('parameters','var') || isempty(parameters)
    parameters.experiment.name='2012-09-28_ex145_07x2_ROI2017';
    parameters.scale.x = '11.24';
    parameters.scale.y = '11.24';
    parameters.scale.z = '28';
    parameters.offset.x = '0';
    parameters.offset.y = '0';
    parameters.offset.z = '0';
end
if ~exist(outputFolder,'dir')
    mkdir(outputFolder);
end
if ~exist('maxSegId','var') || isempty(maxSegId)
    maxSegId = max(cat(1,agglo(:).nodes),[],1);
    maxSegId = maxSegId(4);
end
% Set colors to be used
colors = distinguishable_colors(min(length(agglo),100), [0 0 0; 1 1 1]);
colors = repmat(colors, ceil(length(agglo)/100), 1);
colors(:,4) = 0;
allSkel = skeleton();
allSkel.parameters = parameters;
allSkel.scale = structfun(@str2double, allSkel.parameters.scale)';

for tr=1:length(agglo)
    if ~isempty(agglo(tr).nodes)
        % Generate parameters for skeleton
        theseCoM = agglo(tr).nodes(:,1:3);
        skel = skeleton();
        if ~isfield(agglo(tr),'comments')
            agglo(tr).comments = [];
        end
        skel = skel.addTree([treeNames{tr} '_' num2str(size(theseCoM,1))],[theseCoM,ones(size(agglo(tr).nodes,1),1)*10],agglo(tr).edges,colors(tr,:),[],agglo(tr).comments);
        allSkel = allSkel.addTreeFromSkel(skel);
        if ~exist('filename','var') || isempty(filename)
            skel.parameters = parameters;        
            skel.scale = structfun(@str2double, skel.parameters.scale)';
            skel.write(fullfile(outputFolder,[ treeNames{tr} '.nml']), [], 1);
        end
    end
end
if exist('filename','var') && ~isempty(filename)
    allSkel.write(fullfile(outputFolder,filename), [], 1)
end
mappingFile = fullfile(outputFolder, [treeNames{1} '.txt']);
script = WK.makeMappingScript(maxSegId, Superagglos.transformAggloNewOldRepr(agglo), false);
fileHandle = fopen(mappingFile, 'w');
fwrite(fileHandle, script);
fclose(fileHandle);

end
