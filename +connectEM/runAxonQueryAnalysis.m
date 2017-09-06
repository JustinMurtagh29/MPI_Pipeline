function runAxonQueryAnalysis(state,superagglo,comment,queries,chiasma)

if nargin < 2
    superagglo = false;
end

if nargin < 3
    comment = false;
end

if nargin < 4
    queries = false;
end

if nargin < 5
    chiasma = false;
end

param = load('/gaba/u/mberning/results/pipeline/20170217_ROI/allParameterWithSynapses.mat');
param = param.p;

% Optional after chiasmata splitting: 
if chiasma
    tic;connectEM.generateAxonEndings(param, state);toc;
end

if ~comment
    tic;connectEM.getAggloQueryOverlapA(param, state);toc;
end

if comment
    tic;connectEM.getAggloQueryOverlapB_comment(param, state);toc;
else
    tic;connectEM.getAggloQueryOverlapB(param, state);toc;
end

tic;connectEM.flightEndingOverlapRun(param, state);toc;

if ~comment
    tic;connectEM.makeEndingCaseDistinctions(param, state);toc;
end

if superagglo
    tic;connectEM.createNewSuperagglos(param,state);toc;
end

if queries
    tic;connectEM.generateAxonQueries(param, state);toc;
end