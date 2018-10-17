outDir = '/tmpscratch/mbeining/L4/';
dendState = 'dendrites_03_v2.mat';
axState = 'axons_03.mat';


load('/gaba/u/mberning/results/pipeline/20170217_ROI/allParameterWithSynapses.mat');
p.seg.root = '/gabaghi/wKcubes_archive/2012-09-28_ex145_07x2_ROI2017/segmentation/1/';
p.seg.backend = 'wkwrap';
aggloFolder = fullfile(p.saveFolder, 'aggloState');

folder = '/gaba/u/mbeining/code/pipeline/+connectEM/+WholeCell/+Data/border-cells_axon-dendrites-split';
files = dir(fullfile(folder,'*.nml'));
filenames = arrayfun(@(x) fullfile(folder,x.name),files,'uni',0);
folder = '/gaba/u/mbeining/code/pipeline/+connectEM/+WholeCell/+Data/center-cells_axon-dendrites-split';
files = dir(fullfile(folder,'*.nml'));
filenames = cat(1,filenames,arrayfun(@(x) fullfile(folder,x.name),files,'uni',0));
% load dendrites and axons
load(fullfile(aggloFolder,dendState),'dendrites','indBigDends');
load(fullfile(aggloFolder,axState),'axons','indBigAxons');
axons = axons(indBigAxons);
dendrites = dendrites(indBigDends);
[axonLUT,axonSegIds] = Superagglos.buildLUT(axons);
[dendriteLUT,dendriteSegIds] = Superagglos.buildLUT(dendrites);
load(fullfile(aggloFolder,'wholeCells_07.mat'))
[wcLUT,wcSegIds] = Superagglos.buildLUT(wholeCells);
dendrites = Superagglos.transformAggloNewOldRepr(dendrites);
axons = Superagglos.transformAggloNewOldRepr(axons);


%% make luts the same size
maxSegId = max([max(axonSegIds),max(dendriteSegIds),max(wcSegIds)]);
if numel(axonLUT) < maxSegId
    axonLUT(maxSegId) = 0;
end
if numel(dendriteLUT) < maxSegId
    dendriteLUT(maxSegId) = 0;
end
if numel(wcLUT) < maxSegId
    wcLUT(maxSegId) = 0;
end
%%
numCoords = zeros(numel(filenames),1);
allSkelCoords = zeros(0,3);
for f = 1:numel(filenames)
    skel = skeleton(filenames{f});
    
    skel = skel.deleteTrees(cellfun(@numel,skel.nodes)/4==0); % delete zero node trees
    skelCoords = cell2mat(cellfun(@(x) x(:,1:3),skel.nodes,'uni',0));
    numCoords(f) = size(skelCoords,1);
    allSkelCoords = cat(1,allSkelCoords,skelCoords);  % putting all skel nodes together
end
warning('OFF','auxiliaryMethods:readKnossosCube')
skelSegIds = Seg.Global.getSegIds(p,allSkelCoords);  % extract the seg Ids of all skel nodes
warning('ON','auxiliaryMethods:readKnossosCube')

skelSegIdsPerFile = cellfun(@nonzeros,mat2cell(skelSegIds,numCoords),'uni',0);

%%
agglos = cell(numel(filenames),1);
allOutDir5 = agglos;
allOutDir1 = agglos;
for f = 1:numel(filenames)
    wcInd = mode(wcLUT(skelSegIdsPerFile{f}));
    
    % avoid pollution of isosurfaces by small axonic stuff which was already
    % taken into account by dendrites
    overlappingAgglosDend = unique(nonzeros(dendriteLUT(skelSegIdsPerFile{f})));
    axonLUT(ismember(dendriteLUT,overlappingAgglosDend)) = 0;
    overlappingAgglosAx = unique(nonzeros(axonLUT(skelSegIdsPerFile{f})));
    
    agglos{f} = [dendrites(overlappingAgglosDend);axons(overlappingAgglosAx)];
    [~,ia] = sort(cellfun(@numel,agglos{f}),1,'descend');
    agglos{f} = agglos{f}(ia);
    allOutDir5{f} = fullfile(outDir,sprintf('isoForWC%d_%g',wcInd,0.5));
    if ~exist(allOutDir5{f},'dir')
        mkdir(allOutDir5{f});
    end
    allOutDir1{f} = fullfile(outDir,sprintf('isoForWC%d_%g',wcInd,0.1));
    if ~exist(allOutDir1{f},'dir')
        mkdir(allOutDir1{f});
    end
end

Visualization.exportAggloToAmira(p,agglos,allOutDir5,'reduce',0.5,'smoothSizeHalf',4,'smoothWidth',8);
Visualization.exportAggloToAmira(p,agglos,allOutDir1,'reduce',0.1,'smoothSizeHalf',4,'smoothWidth',8);

clearvars dendrites axons axonLUT dendriteLUT agglos overlappingAgglosDend overlappingAgglos Ax ia skelSegIdsPerFile skelSegIds
info = Util.runInfo;
save(fullfile(outDir,'info.mat'),'info')


