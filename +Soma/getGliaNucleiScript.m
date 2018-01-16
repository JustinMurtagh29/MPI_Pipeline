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
    save(outFile, 'info', 'gliaNuclei')
else
    Util.log('File %s already exists and will not be overwritten.');
end
