% extractCandidateAIS
% script to cut out all known soma (outgrown) from a given dendrite state
% (originally on dendrites_02 (wo myelin fix) and Robin's soma)
% -------------------------------------------------------------------------
% Author: Matej Zecevic <matej.zecevic@brain.mpg.de>

%% loading state & lookup table
% orig state to be used
% ---------------------
tic;
origStatePath = '/home/zecevicm/Desktop/connectomics_git/L4_apicalDendrites/getClassificationOfComponents/dendrites_02.mat';
load(origStatePath);
agglosNew = dendrites;
%clear dendrites;
agglos = connectEM.transformAggloNewOldRepr(dendrites);
metaPath = '/mnt/gaba/gaba/u/mberning/results/pipeline/20170217_ROI/segmentMeta.mat';
sM = load(metaPath);
disp('finished loading state.'); toc;

% get soma to cut out
% -------------------
tic;
% load somas from Robin
load('/mnt/gaba/gaba/u/mberning/results/pipeline/20170217_ROI/aggloState/somas.mat');
somaAgglos = somas(:,3);
clear somas;

% awesome code by Benedikt
% idToSoma gives mapping from segment id to soma (linear index)
% isInSoma cell array with all soma ids intersecting with given agglo
idToSoma(cell2mat(somaAgglos)) = repelem(1:length(somaAgglos), cellfun(@length, somaAgglos));
m = max(cellfun(@max, agglos));
idToSoma(end+1:m) = 0;
isInSoma = cellfun(@(x)nonzeros(unique(idToSoma(x))), agglos, 'uni', 0);

% collect the agglos which actually contain a soma (some contain even two)
h = cellfun(@(x) ~isempty(x), isInSoma, 'uni', 0);
intersectingInd = find(vertcat(h{:}) == 1);
% count multiple occurences, actually only 6 cases are double out of 1331
% when using the dendrites_02.mat
intersectingVal = arrayfun(@(x) isInSoma{x}, intersectingInd, 'un', 0);
intersectingNum = arrayfun(@(x) numel(isInSoma{x}), intersectingInd, 'un', 0);
intersectingNum = vertcat(intersectingNum{:});
intersecting = horzcat(intersectingInd, intersectingNum);
disp('finished collecting soma intersections.'); toc;

% load graph neighbours
tic;
load('/mnt/gaba/gaba/u/mberning/results/pipeline/20170217_ROI/graphNew.mat','neighbours');
disp('finished loading graph neighbours.'); toc;

%check for 204269, 204873,205721 - 89 and 204815,204881,204898 - 80
%% test & write out

for n = 1:size(somaAgglos,1)
    outgrownSoma{n} = unique(cat(1,somaAgglos{n},cat(2,neighbours{somaAgglos{n}})'));
end

agglosToWrite = {agglos{204269},agglos{204873},agglos{205721},somaAgglos{89},...
    unique(cat(1,somaAgglos{n},cat(2,neighbours{somaAgglos{n}})'))};
writePath = '/home/zecevicm/Desktop/connectomics_git/L4_apicalDendrites/20170912_somaNewTry_agglo26.nml';
distthr = 3000;
%Apicals.writeAgglos( agglosToWrite, sM, writePath);                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          ', distthr );
tic;
connectEM.generateSkeletonFromNodes(writePath,...
    cellfun(@(x) sM.point(:,x)', agglosToWrite,'uni',0), ...
    arrayfun(@(x) strcat('skelNod_',num2str(x,'%.2i')), 1:numel(agglosToWrite),'uni',0),[],[],[],distthr); 
fprintf('created %d skeletons from filtered agglos.\n', length(agglosToWrite)); toc;



% mask out the intersecting ones because they will be
% there as components in new state
maskRest = ~vertcat(h{:});


% cut out soma to get all components
% ----------------------------------
tic;
% load segment meta first
metaPath = '/mnt/gaba/gaba/u/mberning/results/pipeline/20170217_ROI/segmentMeta.mat';
sM = load(metaPath);
componentState = {};
for i=1:length(intersectingInd)
    ind = intersectingInd(i);
    % for cases where there are 2 soma (maybe even more)
%     if length(isInSoma{ind}) > 1
%         % cut first and then use cut version for others (not sure this will
%         % work)
%         [cut, mask] = Apicals.cutOutSoma(sM, somaAgglos{isInSoma{ind}(1)}, agglos{ind}, 700);
%         componentState = vertcat(componentState, cut);
%         for j=1:length(isInSoma{ind})
%             [cut, mask] = Apicals.cutOutSoma(sM, somaAgglos{isInSoma{ind}(j)}, cut, 700);
%             componentState = vertcat(componentState, cut);
%         end
%     else
    if ~isempty(isInSoma{ind}) && length(isInSoma{ind}) == 1
        [cut, mask] = Apicals.cutOutSoma(sM, somaAgglos{isInSoma{ind}}, agglos{ind}, 700);
        
        % function to delete nodes
        aggloAdapted = Apicals.deleteNodesFromAgglo(agglosNew(ind), mask);
        % collect trees
        trees = Graph.findConnectedComponents(aggloAdapted.edges); % -> results in soma (need to convert nodes to segids to work for wK)
        % convert node to segid
        aggloComponents = [];
        for t=1:length(trees)
           aggloComponents = vertcat(aggloComponents, {agglos{ind}(trees{t})});
        end
        
        componentState = vertcat(componentState, aggloComponents);
    end
end
disp('-----------------------------------');
disp('finished computing component state.'); toc;

% connect all together and save
% -----------------------------
componentWithRestState = vertcat(componentState, agglos(maskRest));
% save state
% save('/home/zecevicm/Desktop/connectomics_git/L4_apicalDendrites/20170912_componentState_test.mat', 'componentState');

