% script for the soma target class gallery creation
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

info = Util.runInfo();


%% get soma agglos

Util.log('Loading data');
p = Gaba.getSegParameters('ex145_ROI2017');
p.seg.backend = 'wkwrap';
p.seg.root = ['/gabaghi/wKcubes_archive/2012-09-28_ex145_07x2_ROI2017/' ...
    'segmentation/1/'];
m = load(p.agglo.somaFile);
somaAgglos = m.somaAgglos(m.centerSomaIdx, 1);


%% soma synapses

% see script Soma.somaSynapses in the pipeline


%% amira issfs of somata

outputFolder = '/tmpscratch/bstaffle/data/L4/somata/issfs_amira/';
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder)
end
save('RunInfo.mat', 'info');
Util.log('Calculating isosurfaces and storing them at %s.', outputFolder);
Visualization.exportAggloToAmira(p, somaAgglos, outputFolder, ...
    'reduce', 0.2, 'smoothSizeHalf', 4, 'smoothWidth', 8);
Util.log('Finished isosurface calculation.');
