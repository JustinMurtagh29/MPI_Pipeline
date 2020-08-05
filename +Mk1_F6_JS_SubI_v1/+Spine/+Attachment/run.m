% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
% Modified by
%   Sahil Loomba <sahil.loomba@brain.mpg.de>

clear;

%% configuration
rootDir = '/tmpscratch/sahilloo/data/Mk1_F6_JS_SubI_v1/pipelineRun_mr2e_wsmrnet/';
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

nmlDir = fullfile(param.saveFolder, 'tracings', 'spine-attachment');
dendFile = fullfile(param.saveFolder, 'aggloMat', '20191227T134319_ga_20191224T235355optimParams_agglomeration/20191227T220548_results.mat');

aggloParam = struct;
% NOTE(amotta): The spine head probability was optimized using:
% L23.TypeEM.buildSegmentClassifier
% git@gitlab.mpcdf.mpg.de:connectomics/amotta.git 6a46f4500a4392aab1a5c33b8f0b4f496d1f22b5 (dirty)
% amotta@gaba01. MATLAB 9.3.0.713579 (R2017b). 11-Dec-2019 15:19:38
%
% More detailed description in
% L23/2018-10-17-Spine-Head-Detection/notes.md
aggloParam.minSpineHeadProb = 0.661619; %0.5545;
aggloParam.maxVesselScore = 0.5;
aggloParam.maxNucleusScore = 0.5;

% NOTE(amotta): Spine head agglomerates generated with the following edge
% probability threshold were exported to and inspected in webKnossos with inspectSHAgglos.m
aggloParam.minEdgeProb = 0.98;

% NOTE(amotta): See commit 3f82cec80114079534b1ce10e44f9af7ab2af3cc for a
% validation of these parameters. This was documented in:
% L23/2018-10-23-Spine-Attachment/notes.md
% NHP: commit 8da753d0aace68ea7a147bc80246df50fc3abcba
attachParam = struct;
attachParam.minEdgeProb = 0.1062;%0.0041;
attachParam.maxNumSteps = round(15.8081);%round(12.20);
attachParam.maxAxonProb = 1;%0.9510;
attachParam.minDendProb = 0.0015;%0.1659;
attachParam.maxAstroProb = 0.9294;%0.5434;

info = Util.runInfo();

%% complete configuration
[outDir, outFile] = fileparts(dendFile);
outFile = sprintf('%s_auto-spines_v3.mat', outFile);
outFile = fullfile(outDir, outFile);
clear outDir;

Util.log('Reading dendrites from "%s"', dendFile);
Util.log('Writing results to "%s"', outFile);

%% actually do some work
Util.log('Building spine heads');
sh = struct;
[sh.agglos, sh.edges] = ...
    L4.Spine.Head.buildAgglos(param, aggloParam);

Util.log('Loading graph');
graph = fullfile(rootDir, 'graph.mat');
graph = load(graph, 'edges', 'prob');
graph = Graph.addNeighbours(graph);

Util.log('Loading dendrites');
dendIn = load(dendFile);
dendIn.dendrites = newToOldSuperAgglo(dendIn.dendrites);
dendIn.dendrites = SuperAgglo.clean(dendIn.dendrites);

Util.log('Running spine head attachment');
trunkMask = true(numel(dendIn.dendrites), 1);
[sh.neck, sh.attached, out] = L4.Spine.Head.attach( ...
    param, attachParam, graph, sh, dendIn, trunkMask);

% TODO(amotta): Enable check when spine attachment is fixed
out.dendrites = SuperAgglo.clean(out.dendrites, false);

%% Completing and writing output
Util.log('Completing output');
    
out.shAgglos = sh.agglos;
out.shEdges = sh.edges;
out.edges = sh.neck;
out.attached = sh.attached;
out.info = info;

Util.log('Writing output');
Util.saveStruct(outFile, out);
Util.protect(outFile);

Util.log('Done');

%% Utilities
function agglos = newToOldSuperAgglo(agglos)
    for curIdx = 1:numel(agglos)
        agglos(curIdx).nodes = horzcat(agglos(curIdx).nodes, ...
            cast(agglos(curIdx).segIds, 'like', agglos(curIdx).nodes));
    end
    
    agglos = rmfield(agglos, 'segIds');
end
