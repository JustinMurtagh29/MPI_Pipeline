% shDensityCompScript 
% script to get the spine head densities for the components of a given agglo
% by cutting out the corresponding soma first - identify AIS (designed for L4)
% -------------------------------------------------------------------------
% Author: Matej Zecevic <matej.zecevic@brain.mpg.de>

%% load necessary files

dendAggloStatePath = '/gaba/u/mberning/results/pipeline/20170217_ROI/aggloState/dendrites_03.mat';
spineHeadsStatePath = '/tmpscratch/zecevicm/L4/Apicals/20170829_spineheadsAttached.mat';
dendLensPath = '/tmpscratch/zecevicm/L4/Apicals/20170829_DendriteAggloLengths.mat';
parameterPath = '/gaba/u/mberning/results/pipeline/20170217_ROI/allParameterWithSynapses.mat';
metaPath = '/gaba/u/mberning/results/pipeline/20170217_ROI/segmentMeta.mat';
local = '';
[ agglosNew, agglos, agglos_all, points, indBigDends, bounds, spineHeads, dendLens, sM, p ] = ...
    Apicals.loadDendAggloState(dendAggloStatePath, spineHeadsStatePath, dendLensPath, parameterPath, metaPath, local);

load('/tmpscratch/zecevicm/L4/Apicals/20170829_spineheadCounts(OnlyBig).mat', 'shCounts');


%% get densities

% collect the densities
shDensity = shCounts ./ dendLens;
% replace NaN with zero spine density
shDensity(isnan(shDensity)) = 0;

%% plots spine density per um

% %figure;
% h=histogram(shDensity, 30);
% xlabel('Spine density: # spines per um');
% ylabel('# dendrite agglos');
% title(['Spine Head density for ' num2str(length(shDensity)) ' agglos [Bin width: ' num2str(h.BinWidth) ']']);
% %hold on

%% get somas

% load somas from Robin
load('/gaba/u/rhesse/forBenedikt/somasNoClosedHoles.mat');
somaAgglos = somas(:,3);
% ones he probably excluded?
somaAgglos(94) = []; somaAgglos(93) = []; somaAgglos(91) = [];
somaAgglos(90) = []; somaAgglos(86) = []; somaAgglos(75) = [];

%% create lookup table

% awesome code by Benedikt
% idToSoma gives mapping from segment id to soma (linear index)
% isInSoma cell array with all soma ids intersecting with given agglo
idToSoma(cell2mat(somaAgglos)) = repelem(1:length(somaAgglos), cellfun(@length, somaAgglos));
m = max(cellfun(@max, agglos));
idToSoma(end+1:m) = 0;
isInSoma = cellfun(@(x)nonzeros(unique(idToSoma(x))), agglos, 'uni', 0);

%% get the cutoff agglos and look at density distribution

% get one example
%cut = Apicals.cutOutSoma(sM, somaAgglos{11}, agglos{3}, 700);
% 16 is now to 78, which was nice kinda - know what is AIS through Ali
[cut, mask] = Apicals.cutOutSoma(sM, somaAgglos{78}, agglos{16}, 700);

% collect the agglos which actually contain a soma (some contain even two)
h = cellfun(@(x) ~isempty(x), isInSoma, 'uni', 0);
intersectingInd = find(vertcat(h{:}) == 1);

%% write 

% agglosToWrite = vertcat(agglos(16), somaAgglos(78), cut);
% 
% tic;
% distthr = 3000;  % maximum distance between nodes that can be connected
% connectEM.generateSkeletonFromNodes('/home/zecevicm/Desktop/connectomics_git/L4_apicalDendrites/20170830_cutoutNiceOne_test.nml',...
%     cellfun(@(x) sM.point(:,x)', agglosToWrite,'uni',0), ...
%     arrayfun(@(x) strcat('skelNod_',num2str(x,'%.2i')), 1:numel(agglosToWrite),'uni',0),[],[],[],distthr); 
% fprintf('created %d skeletons from filtered agglos.\n', length(agglosToWrite)); toc;


%% collect components

% function to delete nodes
aggloAdapted = Apicals.deleteNodesFromAgglo(agglosNew(16), mask);

% collect trees
trees = Graph.findConnectedComponents(aggloAdapted.edges); % -> results in soma (need to convert nodes to segids to work for wK)

% convert node to segid
aggloComponents = [];
for t=1:length(trees)
   aggloComponents = vertcat(aggloComponents, {agglos{16}(trees{t})});
end

agglosToWrite = aggloComponents;

%% get spines attached to components

% collect the segments reached
[ segmentsReached ] = Apicals.collectSegmentsReached(spineHeads.edges, spineHeads.shAgglos);

% check for agglo all the spine heads that got attached
% sanity check: [ shAttached ] = Apicals.attachedTo( agglos{16}, segmentsReached, 16, shCounts );
aggloCompCounts = zeros(length(aggloComponents),1);
for i=1:length(aggloComponents)
    [ shAttached ] = Apicals.attachedTo( aggloComponents{i}, segmentsReached);
    aggloCompCounts(i) = numel(shAttached);
end
% not exactly because I dont have a perfect soma cut out
% assert(sum(aggloCompCounts) + countForSoma == shCounts(16));


%% get densities for comp

% get lengths for components
% use Alessandros function
lensComp = Agglo.calcPathLengths(p, aggloComponents);
lensComp = lensComp ./ 1E3; % um

% get densities
shdComp = aggloCompCounts ./ lensComp; %per um
% replace NaN with zero spine density
shdComp(isnan(shdComp))=0; %zero spines

% get segment counts
l = cellfun(@length, aggloComponents, 'un', 0);
l = vertcat(l{:});

% create table
compTable = table([1:length(aggloComponents)]',l,lensComp,aggloCompCounts,shdComp,...
    'VariableNames',{'TreeNumber' 'SegmentsCount' 'Pathlength' 'SpineHeadsCount' 'SHdensity'});

%% plot

% %figure;
% histogram(shdComp, 60);
% xlabel('Spine density: # spines per um');
% ylabel('# dendrite agglos');
% title(['Spine Head density for ' num2str(length(aggloComponents)) ' agglos filtered']);
% %hold on

