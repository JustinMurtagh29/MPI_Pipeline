aggloName = 'dendrites_wholeCells_03_v2.mat';
mainOutDir = ['/tmpscratch/mbeining/isos/',connName(1:end-4)];
isoParam = {'smoothWidth',5,'smoothSizeHalf',5,'reduce',0.1};  % heiko suggests 9 , 9 , 0.2

%% load stuff
load('/gaba/u/mberning/results/pipeline/20170217_ROI/allParameterWithSynapses.mat');
p.seg = struct;
p.seg.root = '/tmpscratch/amotta/l4/2012-09-28_ex145_07x2_ROI2017/segmentation/1';
p.seg.backend = 'wkwrap';
disp('Parameters loaded');
outputFolder = fullfile(p.saveFolder, 'aggloState');

if ~exist('graph','var') || ~all(isfield(graph,{'edges','borderIdx'}))
    graph = load(fullfile(p.saveFolder, 'graphNew.mat'),'edges','borderIdx');
end

load(fullfile(p.saveFolder,'aggloState',aggloName))


disp('all data loaded')
%% export WCs iso
Visualization.exportAggloToAmira(p, Superagglos.transformAggloNewOldRepr(dendrites(indWholeCells)), mainOutDir, isoParam{:})

