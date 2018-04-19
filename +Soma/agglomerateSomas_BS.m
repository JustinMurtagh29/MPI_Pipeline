%% script to run soma agglomeration
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>
% Based on code by Robin Hesse

info = Util.runInfo(false);
doNucleiDetection = false; %#ok<*UNRCH>
ver = '06';

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
    MouseRoi2016.detectNucleiinBoxes2();

    % bboxes for nuclei (bugfix for wrong bbox in agglomeration)
    p = Gaba.getSegParameters('ex145_ROI2017');
    m = load(fullfile(p.saveFolder, 'soma', 'NucleiCoordinates.mat'), 'rp');
    rp = m.rp;
    bboxes = Soma.nucleiBbox(rp, p);
    nucMaskFolder = '/gaba/u/mberning/results/pipeline/20170217_ROI/soma/Nuclei/';
    for i = 1:length(rp)
        filepath = fullfile(nucMaskFolder, sprintf('Nucleus%d.mat', i));
        m = load(filepath);
        if all((diff(bboxes{i}, [], 2) + 1) == size(m.nucleus)')
            m.bbox = bboxes{i};
            save(filepath, '-struct', 'm');
            clear m
        else
            warning('Nucleus %d bbox does not the the mask.', i);
        end
    end

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
% somaLIst = fullfile('+Soma/somaDoc_BS.xlsx');
numNuclei = 125;

% choose somaIds
somaIDs = 1:numNuclei;

% get glia
[~, ~, gliaIds] = xlsread(somaList, sprintf('B2:B%d', numNuclei + 1));
gliaIds = find(~cellfun(@isnan, gliaIds));

% get somata in center
[~, ~, centerSomaIdx] = xlsread(somaList, sprintf('C2:C%d', numNuclei + 1));
centerSomaIdx = ~cellfun(@isnan, centerSomaIdx);

%remove glia
somaIDs(gliaIds) = [];
centerSomaIdx(gliaIds) = [];


%% discard edges into myelin/vessel

m = load(p.svg.heuristicFile, 'vesselScore', 'myelinScore');
excludeId = find(m.vesselScore > 0.5 | m.myelinScore > 0.5);
graph.prob(any(ismember(graph.edges, excludeId), 2)) = 0;


%% soma agglo runs

Util.log('Starting soma agglomeration.');
somaAgglos = Soma.getSomaNodesPar(p, graph, meta, rp, gb, 0.98, 2000, ...
    somaIDs);
somaAgglos = somaAgglos(:,[3 2]);
% somaAgglos2 = Soma.getSomaNodesPar(p, graph, meta, rp, gb, 0.98, 2000, ...
%     somaIDs);
% somaAgglos3 = Soma.getSomaNodesPar(p, graph, meta, rp, gb, 0.99, 2000, ...
%     somaIDs);


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


%% add segments that are mostly surrounded by agglos

fractT = 0.8;
iter = 10;
toAddSingleSegs = Soma.addSurroundedSegments(somaAgglos(:,1), ...
    graph.edges(~isnan(graph.borderIdx), :), gb.borderSize, fractT, ...
    graph.edges(isnan(graph.borderIdx), :), iter);

somaAgglos_woSurSegs = somaAgglos;
for i = 1:size(somaAgglos, 1)
    somaAgglos{i,1} = cat(1, somaAgglos{i,1}(:), toAddSingleSegs{i}(:));
end

%% save the automated results

outFile = fullfile(outputFolder, sprintf('somaAgglo_%s.mat', ver));
if ~exist(outFile, 'file')
    Util.log('Storing automated agglomeration results at %s.', outFile)
    save(outFile, 'somaAgglos', 'somaAgglos_woSurSegs', ...
        'centerSomaIdx', 'info')
else
    Util.log('File %s already exists and will not be overwritten.');
end


%% convert to superagglos

p = Gaba.getSegParameters('ex145_ROI2017');
outFile = fullfile(p.agglo.saveFolder, ['somata_', ver, '.mat']);
assert(~exist(outFile, 'file'));
Soma.agglosToSuperagglos(p);

