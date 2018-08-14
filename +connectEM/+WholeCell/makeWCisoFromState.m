function makeWCisoFromState(state,reduce,outDir)

if ~exist('reduce','var')
    reduce = 0.1;
end
if ~exist('outDir','var')
    outDir = '/tmpscratch/mbeining/L4/';
end

load('/gaba/u/mberning/results/pipeline/20170217_ROI/allParameterWithSynapses.mat');
outputFolder = fullfile(p.saveFolder, 'aggloState');
m = load(fullfile(outputFolder,[state,'.mat']));

if ~isfield(m,'WholeCellId')
    %% load soma whole cell agglos
    somaAgglos = connectEM.getSomaAgglos(fullfile(outputFolder,'somas_with_merged_somas.mat'),'center');
    somaAgglos = somaAgglos(~cellfun(@isempty,somaAgglos)); % remove not found somata
    somaSegIds = cell2mat(somaAgglos);
    % remove duplicate segIds
    [~,ic] = unique(somaSegIds);
    duplicates = somaSegIds(setdiff(1:numel(somaSegIds),ic));
    somaAgglos = cellfun(@(x) setdiff(x,duplicates),somaAgglos,'uni',0);
    somaSegIds = cell2mat(somaAgglos);
    somaLUT(somaSegIds) = repelem(1:numel(somaAgglos),cellfun(@numel,somaAgglos));
    disp('Soma whole cell agglos loaded')
    
    [dendriteLUT,dendriteSegIds] = Superagglos.buildLUT(dendrites);
    [ismem,ind] = ismember(somaSegIds,dendriteSegIds);
    % get each dend id which contains most of the seg ids of each soma
    m.WholeCellId = accumarray(somaLUT(somaSegIds(ismem))',dendriteLUT(dendriteSegIds(ind(ismem)))',[],@mode);
end

outDir = fullfile(outDir,sprintf('isoForWC_%s_%g',state,reduce));

if ~exist(outDir,'dir')
    mkdir(outDir);
end

Util.log('Loading stuff..')
Util.log('Now running main viz')
p.seg.root = ['/gabaghi/wKcubes_archive/2012-09-28_ex145_07x2_ROI2017/segmentation/1/'];
p.seg.backend = 'wkwrap';
Visualization.exportAggloToAmira(p,m.dendrites(m.WholeCellId),outDir,'reduce',reduce,'smoothSizeHalf',4,'smoothWidth',8);
info=Util.runInfo(false);
save(fullfile(outDir,'info.mat'),'info');
