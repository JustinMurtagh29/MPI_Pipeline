function makeIso(agglos,outDir,p)
Util.log('initialize stuff')
if ~exist(outDir,'dir')
    mkdir(outDir);
end
if ~exist('p','var')
    load('/gaba/u/mberning/results/pipeline/20170217_ROI/allParameterWithSynapses.mat');
    p.seg.root = ['/gabaghi/wKcubes_archive/2012-09-28_ex145_07x2_ROI2017/segmentation/1/'];
    p.seg.backend = 'wkwrap';
end
Util.log('Now running main viz')
Visualization.exportAggloToAmira(p,agglos,outDir,'reduce',0.5,'smoothSizeHalf',4,'smoothWidth',8);
info=Util.runInfo(false);
save(fullfile(outDir,'info.mat'),'info','agglos');