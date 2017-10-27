%% script to run soma agglomeration
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>
% Based on code by Robin Hesse

info = Util.runInfo(false);
doNucleiDetection = false; %#ok<*UNRCH>
ver = '03';

outputFolder = sprintf(['/gaba/u/mberning/results/pipeline/' ...
    '20170217_ROI/soma_BS/']);

if ~exist(outputFolder, 'dir')
    mkdir(outputfolder)
end

% cross check if file already exists
if exist(fullfile(outputFolder, sprintf('cells_%s.mat', ver)), 'file')
    warning(['Found files associated with version %s. Results of this ' ...
        'run will not be stored automatically.'], ver);
end

%% nuclei detection (just as a reference; wasn't run like this)
% requires the repo nuclear_pores on the path

if doNucleiDetection

    % nuclei coordinates and bbox
    rp = NuclearPores.mag4NucleiDetection(); 
    outFile = ['/gaba/u/mberning/results/pipeline/20170217_ROI/soma/' ...
        'NucleiCoordinates.mat'];
    if ~exist(outFile, 'file')
        save(outFile,'rp');
    end
    
    % nuclei mag1 masks (output files need to be copied to
    % /gaba/u/mberning/results/pipeline/20170217_ROI/soma/Nuclei/
    MouseROI2017.detectNucleiinBoxes2();
    
end


%% load svg data

Util.log('Loading SVG data.');
p = Gaba.getSegParameters('ex145_ROI2017');
graph = load([p.saveFolder 'graphNew.mat'], 'prob', 'edges', 'borderIdx');
meta = load(p.svg.segmentMetaFile, 'segIds', 'point', 'voxelCount');
m = load(fullfile(p.saveFolder, 'soma', 'NucleiCoordinates.mat'), 'rp');
rp = m.rp;
gb = load(p.svg.borderMetaFile, 'borderCoM', 'borderSize');


%% select somata

thisFolder = fileparts(mfilename('fullpath'));
somaList = fullfile(thisFolder, 'somaDoc_BS.xlsx');
numNuclei = 125;

% choose somaIds
somaIDs = 1:numNuclei;

% get glia
[~, ~, gliaIds] = xlsread(somaList, sprintf('B2:B%d', numNuclei + 1));
gliaIds = find(~cellfun(@isnan, gliaIds));

%remove glia
somaIDs(gliaIds) = [];


%% discard edges into myelin/vessel

m = load(p.svg.heuristicFile, 'vesselScore', 'myelinScore');
excludeId = find(m.vesselScore > 0.5 | m.myelinScore > 0.5);
graph.prob(any(ismember(graph.edges, excludeId), 2)) = 0;


%% soma agglo runs

Util.log('Starting soma agglomeration.');
somaAgglos = Soma.getSomaNodesPar(p, graph, meta, rp, gb, 0.98, 1500, somaIDs);
somaAgglos = somaAgglos(:,[3 2]);


%% remove blood vessels & myelin
% % (should maybe excluded in agglo procedure directly)
% 
% m = load(p.svg.heuristicFile, 'vesselScore', 'myelinScore');
% isVessel = m.vesselScore > 0.5;
% isMyelin = m.myelinScore > 0.5;
% toExclude = isVessel | isMyelin;
% idx = cellfun(@(x)any(isnan(x(:))), somaAgglos(:,1));
% somaAgglos(~idx,1) = cellfun(@(x)x(~toExclude(x)), somaAgglos(~idx,1), ...
%     'uni', 0);


%% save the automated results

outFile = fullfile(outputFolder, sprintf('cells_%s.mat', ver));
if ~exist(outFile, 'file')
    Util.log('Storing automated agglomeration results at %s.', outFile)
    save(outFile, 'somaAgglos', 'info')
else
    Util.log('File %s already exists and will not be overwritten.');
end
