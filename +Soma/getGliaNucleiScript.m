% script to get the glia nuclei (based on agglomerateSomas_BS.m)
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

info = Util.runInfo();
outputFolder = sprintf(['/gaba/u/mberning/results/pipeline/' ...
    '20170217_ROI/soma_BS/']);


%% glia ids

thisFolder = fileparts(mfilename('fullpath'));
somaList = fullfile(thisFolder, 'somaDoc_BS.xlsx');

% get glia
numNuclei = 125;
[~, ~, gliaIds] = xlsread(somaList, sprintf('B2:B%d', numNuclei + 1));
gliaIds = find(~cellfun(@isnan, gliaIds));


%% get glia nuclei

gliaNuclei = Soma.getNuclei(gliaIds);


%% save results

outFile = fullfile(outputFolder, 'gliaNuclei.mat');
if ~exist(outFile, 'file')
    Util.log('Storing automated agglomeration results at %s.', outFile)
    save(outFile, 'info', 'gliaNuclei', 'gliaIds')
else
    Util.log('File %s already exists and will not be overwritten.');
end


%% calculate isosurfaces

p = Gaba.getSegParameters('ex145_ROI2017');
outDirIssfs = '/tmpscratch/bstaffle/data/L4/somata/nuclei_issfs/';
if ~exist(outDirIssfs, 'dir')
    mkdir(outDirIssfs)
end
m = load(outFile);
gliaNuclei = m.gliaNuclei;
Util.log('Calculating isosurfaces and storing them at %s.', outDirIssfs);
Visualization.exportAggloToAmira(p, gliaNuclei, outDirIssfs, ...
    'reduce', 0.2, 'smoothSizeHalf', 4, 'smoothWidth', 8);
Util.log('Finished isosurface calculation.');