function makeIso(state,reduce,outDir,type)

if ~exist('reduce','var')
    reduce = 0.1;
end
if ~exist('outDir','var')
    outDir = '/tmpscratch/mbeining/L4/';
end
if ~exist('type','var')
    type = 'all'; % can be border, center or all
end
Util.log('Loading parameters and state')
load('/gaba/u/mberning/results/pipeline/20170217_ROI/allParameterWithSynapses.mat');
outputFolder = fullfile(p.saveFolder, 'aggloState');
m = load(fullfile(outputFolder,[state,'.mat']));

if ~isfield(m,'WholeCellId')
    Util.log('Overlap state with soma agglos')
    %% load soma whole cell agglos
    %% load all soma whole cell agglos
    somaAgglos = connectEM.getSomaAgglos(fullfile(outputFolder,'somas_with_merged_somas.mat'),type);
    somaAgglos = somaAgglos(~cellfun(@isempty,somaAgglos)); % remove not found somata
    somaSegIds = cell2mat(somaAgglos);
    % remove duplicate segIds
    [~,ic] = unique(somaSegIds);
    duplicates = somaSegIds(setdiff(1:numel(somaSegIds),ic));
    % somaDuplicateIds = cellfun(@(x) any(intersect(x,duplicates)),somaAgglos);
    somaAgglos = cellfun(@(x) setdiff(x,duplicates),somaAgglos,'uni',0);
    somaSegIds = cell2mat(somaAgglos);
    somaLUT(somaSegIds) = repelem(1:numel(somaAgglos),cellfun(@numel,somaAgglos));
    disp('Soma whole cell agglos loaded')
    
    [dendriteLUT,dendriteSegIds] = Superagglos.buildLUT(m.dendrites);
    [ismem,ind] = ismember(somaSegIds,dendriteSegIds);
    % get each dend id which contains most of the seg ids of each soma
    m.WholeCellId = accumarray(somaLUT(somaSegIds(ismem))',dendriteLUT(dendriteSegIds(ind(ismem)))',[],@mode);
end

outDir = fullfile(outDir,sprintf('isoForWC_%s_%s_%g',type,state,reduce));

if ~exist(outDir,'dir')
    mkdir(outDir);
end

Util.log('Now running iso creation')
p.seg.root = ['/gabaghi/wKcubes_archive/2012-09-28_ex145_07x2_ROI2017/segmentation/1/'];
p.seg.backend = 'wkwrap';
Visualization.exportAggloToAmira(p,Superagglos.transformAggloNewOldRepr(m.dendrites(m.WholeCellId)),outDir,'reduce',reduce,'smoothSizeHalf',4,'smoothWidth',8);
info=Util.runInfo(false);
save(fullfile(outDir,'info.mat'),'info');
Util.log('Done')