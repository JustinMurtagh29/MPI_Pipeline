function multipleMappingsForOxalis(p, upperT)

% Find maximum ID in dataset and display (should be added automatically to layer.json)
load([p.local(end,end,end).saveFolder 'localToGlobalSegId.mat']);
numEl = max(globalIds);
display(['Maximum ID in (globalized) segmentation: ' num2str(numEl)]);
    
% Load global correspondences (e.g. how to continue over cube borders, variable components)
load([p.saveFolder 'globalCorrespondences.mat']);

% Create structure & write JSON for correspondence mapping
coMapping.name = 'correspondences';
coMapping.classes = components; % from globalCorrespondences.mat, see findGlobalCC.m for calculation
writeJson([p.saveFolder 'correspondences.json'],coMapping); 

% Iterate over cubes, load all edges and probabilities, apply upper & lower treshold, 
allComponents = {};
for i=1:size(p.local,1)
    for j=1:size(p.local,2)
        for k=1:size(p.local,3)
            load(p.local(i,j,k).edgeFile);
            load(p.local(i,j,k).probFile);
            % Fix at some point: Definition of upper threshold should be done in setParameterSettingsBig.m, done here for easier modification
            joinEdges = prob > upperT;
            % Calculate, save and collect all GP based components for each cube
            edgesToJoin = edges(joinEdges,:);
            componentsNew = findConnectedComponents(edgesToJoin);
            allComponents = [allComponents; componentsNew]; % this line is probably slowing this down substaintially, optimize if problem 
            Util.save([p.local(i,j,k).saveFolder 'joinedComponents.mat'], componentsNew);
        end
    end
end
% Save all GP probabilities
Util.save([p.saveFolder 'joinedComponents.mat'], allComponents);

% Create structure & write JSON for mapping based on Gaussian process probabilities
gpMapping.name = 'gaussianProcess';
gpMapping.classes = allComponents;
writeJson([p.saveFolder 'gaussianProcess.json'], gpMapping); 

end

