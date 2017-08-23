function [ synIdx ] = getSynapsePredictions( synScores, synT, ...
    doDirectionExcl, prob, probT, edges, myelinScore, myelinT )
%GETSYNAPSEPREDICTIONS Get the synapse prediction using exlusions based on
%prob and myelin scores.
% INPUT synScores: [Nx2] float
%           Synapse prediction for an edge.
%       synT: double
%           Synapse score threshold.
%       doDirectionExcl: (Optional) flag
%           Exclude synapses if both directions have a score above -2.
%           (Default: false)
%       prob: (Optional) [Nx1] float
%           Merge continuity probability for each edge.
%           (Default: no exclusion based on prob).
%       probT: (Optional) float
%           Threshold above which synapse preditions are discarded. This
%           can get rid of intracellular synapse FPs (within large
%           processes/somata/at mito borders).
%           (Default: 0.7 if prob is provided)
%       edges: (Optional) [Nx2] int
%           Global edges list.
%           (Only needed if myeline score is used)
%       myelineScores: [Nx1] float
%           Myelin scores for segments.
%           (Default: no exclusion based on myelin).
%       myelinT: float
%           Threshold above which a segment is considered to be at a myelin
%           sheath. Synapses containing segments with myelin are discarded.
%           (Default: 0.375 if myelinScore is provided).
% OUTPUT synIdx: [Nx1] logical
%           Logical indices of predicted synapses.
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

synIdx = any(synScores > synT, 2);

if exist('doDirectionExcl', 'var') && doDirectionExcl
    synIdx = synIdx & ~all(synScores > -2, 2);
end

if exist('prob', 'var') && ~isempty(prob)
    if ~exist('probT', 'var') || isempty(probT)
        probT = 0.70;
    end
    synIdx = synIdx & prob <= probT;
end

if exist('edges', 'var') && exist('myelinScore', 'var') && ...
        ~isempty(edges) && ~isempty(myelinScore)
    if ~exist('myelinT', 'var') || isempty(myelinT)
        myelinT = 0.375;
    end
    synIdx = synIdx & ~any(myelinScore(edges) > myelinT, 2);
end


end

