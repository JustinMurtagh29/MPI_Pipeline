% orig state to be used
% ---------------------
origStatePath = '/home/zecevicm/Desktop/connectomics_git/L4_apicalDendrites/getClassificationOfComponents/dendrites_02.mat';
load(origStatePath);
agglos = Superagglos.transformAggloNewOldRepr(dendrites);

% get soma to cut out
% -------------------
% load somas from Robin
load('/mnt/gaba/gaba/u/rhesse/forBenedikt/somasNoClosedHoles.mat');
somaAgglos = somas(:,3);
% ones he probably excluded?
% get this out at some point probably
somaAgglos(94) = []; somaAgglos(93) = []; somaAgglos(91) = [];
somaAgglos(90) = []; somaAgglos(86) = []; somaAgglos(75) = [];

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

% cut out soma to get all components
% ----------------------------------
% load segment meta first
metaPath = '/mnt/gaba/gaba/u/mberning/results/pipeline/20170217_ROI/segmentMeta.mat';
sM = load(metaPath);
componentState = {};
for i=1:length(intersectingInd)
    ind = intersectingInd(i);
    % for cases where there are 2 soma (maybe even more)
    if length(isInSoma{ind}) > 1
        % cut first and then use cut version for others (not sure this will
        % work)
        [cut, mask] = Apicals.cutOutSoma(sM, somaAgglos{isInSoma{ind}(1)}, agglos{ind}, 700);
        componentState = vertcat(componentState, cut);
        for j=1:length(isInSoma{ind})
            [cut, mask] = Apicals.cutOutSoma(sM, somaAgglos{isInSoma{ind}(j)}, cut, 700);
            componentState = vertcat(componentState, cut);
        end
    elseif ~isempty(isInSoma{ind})
        [cut, mask] = Apicals.cutOutSoma(sM, somaAgglos{isInSoma{ind}}, agglos{ind}, 700);
        componentState = vertcat(componentState, cut);
    end
end

% save state
save('/home/zecevicm/Desktop/connectomics_git/L4_apicalDendrites/20170912_componenState_test.mat', 'componentState');

%% write 

agglosToWrite = componentState; %vertcat(agglos(aggloInd), somaAgglos(78), cut);

tic;
distthr = 3000;  % maximum distance between nodes that can be connected
connectEM.generateSkeletonFromNodes('/home/zecevicm/Desktop/connectomics_git/L4_apicalDendrites/20170912_componenState_test.nml',...
    cellfun(@(x) sM.point(:,x)', agglosToWrite,'uni',0), ...
    arrayfun(@(x) strcat('skelNod_',num2str(x,'%.2i')), 1:numel(agglosToWrite),'uni',0),[],[],[],distthr); 
fprintf('created %d skeletons from filtered agglos.\n', length(agglosToWrite)); toc;
