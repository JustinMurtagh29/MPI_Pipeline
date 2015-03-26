function multipleMappingsForOxalis(p)

% Implement better at some point
load([p.saveFolder 'numEl.mat']);
load([p.saveFolder 'globalCorrespondences.mat']);

% Find maximum ID in dataset
load([p.local(end,end,end).saveFolder 'localToGlobalSegId.mat']);
numEl = max(globalIds);
display(['Maximum ID in (globalized) segmentation: ' num2str(numEl)]);

% Create structure & write JSON for correspondence mapping
coMapping.name = 'correspondences';
coMapping.classes = components; % from globalCorrespondences.mat, see findGlobalCC.m for calculation
writeJson([p.saveFolder 'correspondences.json'], coMapping); 

% Load all edges and probabilities
allComponents = {};
for i=1:size(p.local,1)
    for j=1:size(p.local,2)
        for k=1:size(p.local,3)
            load(p.local(i,j,k).edgeFile);
            load(p.local(i,j,k).probFile);
            % Fix at some point: Definition of upper threshold should be done in setParameterSettingsBig.m
            upperThres = .8;
            joinEdges = prob > upperThres;
            edgesToJoin = edges(joinEdges,:);
            uniqueEdgesToJoin = unique(edgesToJoin);
            componentsNew = findConnectedComponents(edgesToJoin, uniqueEdgesToJoin);
            save([p.local(i,j,k).saveFolder 'joinedComponents.mat'], 'componentsNew');
            allComponents = [allComponents newComponents];
        end
    end
end

% Create structure & write JSON for mapping based on Gaussian process probabilities
gpMapping.name = 'gaussianProcess';
gpMapping.classes = allComponents;
writeJson([p.saveFolder 'gaussianProcess.json'], gpMapping); 

end
