function isSynapse = synScoresToSynEdges(graph, synScore)
    % Use synapse scores and convert to (directed) synapse edge decisions
    % This involves applying a threshold (0) to SynEM scores and
    % applying a threshold (0.5) to neurite continuity scores
    % see: https://mhlablog.net/2016/05/06/getting-rid-of-fp-synapse-detections-using-gp-probabilities/
    % isSynapse is for each edge specified by synScore.edgeIdx

    thisProb = graph.prob(synScore.edgeIdx,:);
    isSynapse = zeros(size(synScore.edgeIdx), 'int8');
    idx = thisProb < 0.5 & (synScore.synScores(:,1) > -1.23 | synScore.synScores(:,2) > -1.23);
    direction = int8((single(synScore.synScores(:,1) > synScore.synScores(:,2)) - 0.5) * 2);
    isSynapse(idx) = direction(idx);
    isSynapse(synScore.synScores(:,2) > 0 & thisProb < 0.5) = -1;

end

