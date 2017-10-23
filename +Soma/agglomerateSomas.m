%% script to run soma agglomeration
% Author: Robin Hesse
% Modified by: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

info = Util.runInfo(false);
doNucleiDetection = false;
ver = '02';

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
load(fullfile('/gaba/u/mberning/results/pipeline/20170217_ROI', 'allParameter.mat'));
graph = load([p.saveFolder 'graphNew.mat'], 'prob', 'edges', 'borderIdx');
meta = load(fullfile(p.saveFolder, 'segmentMeta.mat'), 'segIds', 'point', 'voxelCount');
m = load('/gaba/u/mberning/results/pipeline/20170217_ROI/soma/NucleiCoordinates.mat','rp');
rp = m.rp;
gb = load([p.saveFolder 'globalBorder.mat'], 'borderCoM', 'borderSize');


%% select somata

% choose somaIds
somaIDs = linspace(1,125,125);
%remove glia etc.
somaIDs([6,8,19,21,23,27,28,36,38,41,44,51,57,61,82,87,90,92,95,98,100,104,108,110,117,124]) = [];
somaIDs(61) = []; %makes problems %is actually not even a neuron id75

%% soma agglo runs

Util.log('Starting soma agglomeration.');
cells1 = Soma.getSomaNodesPar(p, graph, meta, rp, gb, 0.98, 1500, somaIDs);

%second run
somaIDs = [5,29,33,34,39,53,55,65,83,85,97,99,101,113,119];
cells2 = Soma.getSomaNodesPar(p, graph, meta, rp, gb, 0.99, 2000, somaIDs);

%third run
somaIDs = [29,33,39,55,83,85,97,99,101,113,119];
cells3 = Soma.getSomaNodesPar(p, graph, meta, rp, gb, 0.995, 2250, somaIDs);

%% save the automated results

outFile = fullfile(outputFolder, sprintf('cells_%s.mat', ver));
if ~exist(outFile, 'file')
    Util.log('Storing automated agglomeration results at %s.', outFile)
    save(outFile, 'cells1', 'cells2', 'cells3', 'info')
else
    Util.log('File %s already exists and will not be overwritten.');
end

%% manual & heuristic corrections

%take corrections from run 2&3
somas = cells1;
somas(5,:) = cells2(1,:);
somas(22,:) = [];
somas(25,:) = cells3(2,:);
somas(26,:) = cells2(4,:);
somas(29,:) = [];
somas(39,:) = cells2(6,:);
somas(41,:) = [];
somas(48,:) = cells2(8,:);
somas(64,:) = cells3(5,:);
somas(66,:) = [];
somas(73,:) = cells3(7,:);
somas(74,:) = cells3(8,:);
somas(75,:) = [];
somas(83,:) = [];
somas(87,:) = [];

%solve special soma
%merges into blood vessel
load('/gaba/u/mberning/results/pipeline/20170217_ROI/heuristicResult.mat');
[ nodes, name, segIds ] = Soma.getSomaNodes( p, graph, meta, rp, gb, 0.98, 1500, 101 );
nodesNew = setdiff(nodes, meta.point(:,find(vesselScore))', 'rows');
segIdsNew = setdiff(segIds, find(vesselScore));
somas{92,1} = nodesNew;
somas{92,2} = name;
somas{92,3} = segIdsNew;


%was classified afterwards as neuron but would caus index shift if added properly
[ nodes, name, segIds ] = Soma.getSomaNodes( p, graph, meta, rp, gb, 0.98, 1500, 82 );
somas{93,1} = nodes;
somas{93,2} = name;
somas{93,3} = segIds;
[ nodes, name, segIds ] = Soma.getSomaNodes( p, graph, meta, rp, gb, 0.99, 2000, 124 );
somas{94,1} = nodes;
somas{94,2} = name;
somas{94,3} = segIds;

outFile = fullfile(outputFolder, sprintf('somas_%s.mat', ver));
if ~exist(outFile, 'file')
    Util.log('Storing automated agglomeration results at %s.', outFile)
    save(outFile, 'somas', 'info');
else
    Util.log('File %s already exists and will not be overwritten.');
end


%add somas that merge into another soma
somas(95,:) = cells1(22,:);
somas(96,:) = cells1(43,:);
somas(97,:) = cells1(69,:);
somas(98,:) = cells1(88,:);
somas(99,:) = cells1(93,:);

outFile = fullfile(outputFolder, ...
    sprintf('somas_with_merged_somas_%s.mat', ver));
if ~exist(outFile, 'file')
    Util.log('Storing automated agglomeration results at %s.', outFile)
    save(outFile, 'somas', 'info');
else
    Util.log('File %s already exists and will not be overwritten.');
end