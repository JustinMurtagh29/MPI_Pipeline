function findGlobalCC(p)

% Collect all global correspondences
files = dir([p.correspondence.saveFolder, '*global.mat']);
correspondences = [];
for i=1:length(files)
    load([p.correspondence.saveFolder files(i).name]);
    correspondences = [correspondences; corrGlobal1 corrGlobal2];
end

% Find connected components
% uniqueIds = unique(correspondences);
[equivalenceClasses, objectClassLabels] = Graph.findConnectedComponents(correspondences);
Util.save([p.saveFolder 'globalCorrespondences.mat'], equivalenceClasses, objectClassLabels);

end

