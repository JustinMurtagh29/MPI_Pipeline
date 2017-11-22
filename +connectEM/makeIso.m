function makeIso(agglos,outDir,reduce)

Util.log('initialize stuff')
if ~exist('reduce','var')
    reduce = 0.2;
end
if ~exist(outDir,'dir')
    mkdir(outDir);
end
if isfield(agglos,'nodes')
    agglos = Superagglos.transformAggloNewOldRepr(agglos);
end
% if ~exist('p','var') || isempty(p)
    load('/gaba/u/mberning/results/pipeline/20170217_ROI/allParameterWithSynapses.mat');
    p.seg.root = '/gabaghi/wKcubes_archive/2012-09-28_ex145_07x2_ROI2017/segmentation/1/';
    p.seg.backend = 'wkwrap';
% end
Util.log('Now running main viz')
Visualization.exportAggloToAmira(p,agglos,outDir,'reduce',reduce,'smoothSizeHalf',4,'smoothWidth',8);
Util.log('PLYs complete. Finished')
info=Util.runInfo(false);
save(fullfile(outDir,'info.mat'),'info','agglos');