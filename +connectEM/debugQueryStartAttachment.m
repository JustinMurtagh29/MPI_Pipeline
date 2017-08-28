% Load all queries from Christians batch that did not attach to start
% agglomerate into variable ff_notAttachingToStart
load('/tmpscratch/mberning/queriesNotAttachingAtStartPosition.mat');
startNodes = cat(1, ff_notAttachingToStart.startNode{:});

% Load all position and directions passed to webKnossos in Christian's
% batch for comparison
% NOTE: Position of queries  on hard disk has changed:
% (see connectEM.visualize...)
batchFolder = '/tmpscratch/scchr/AxonEndings/axonQueryGeneration/EndingsV2/';
for i=1:100
    batch(i) = load(fullfile(batchFolder, ['batch' num2str(i, '%.4i') '.mat']));
end
% Extract from batches again for convenience later (I guess)
q = arrayfun(@(x)x.q, batch, 'uni', 0);
pos = cellfun(@(x)x.pos', q, 'uni', 0);
pos = cat(1, pos{:});
dir = cellfun(@(x)x.dir', q, 'uni', 0);
dir = cat(1, dir{:});
angles = cellfun(@(x)x.pos', q, 'uni', 0);
angles = cat(1, angles{:});
theseAxons = arrayfun(@(x)x.theseAxons, batch, 'uni', 0);
theseAxons = cat(1, theseAxons{:});
clear q i batch batchFolder;

% Load parameter struct of pipeline run
load('/gaba/u/mberning/results/pipeline/20170217_ROI/allParameterWithSynapses.mat');

%% Now rerun the query repositiong script to see what is going wrong
% but only at positions that did not attach (due to start pos) for now
for i=1:length(pos)
    for j=1:length(pos{i})
        % [1 1 1] as startNodes in matlab and pos in webKnossos format
        if any(ismember(startNodes, pos{i}{j} + [1 1 1], 'rows'))
            display(['Query index: ' num2str(i) ', ' num2str(j)]);
            % Calculate node evidence for this start position
            idxInStartNodes = ismember(startNodes, pos{i}{j} + [1 1 1], 'rows');
            idxInNodes = ismember(ff_notAttachingToStart.nodes{idxInStartNodes}, startNodes(idxInStartNodes,:), 'rows');
            segIds = cat(2, ff_notAttachingToStart.segIds{idxInStartNodes}(idxInNodes,:), ...
                ff_notAttachingToStart.neighbours{idxInStartNodes}(idxInNodes,:));
            nodeEvidence = sum(ismember(segIds, theseAxons{i}));
            display(['Node evidence: ' num2str(nodeEvidence)]);
            uniqueSegIds = unique(segIds);
            counts = histc(segIds, uniqueSegIds);
            isInAgglo = ismember(uniqueSegIds, theseAxons{i});
            display(num2str(cat(2, uniqueSegIds', counts', isInAgglo')));
            if nodeEvidence > 13
                display('Skipping. Would have attachted without inner sphere cutting ...');
            else
                % [1 1 1] as pos in webKnossos format but readKnossosRoi
                % expects matlab format
                connectEM.correctQueryLocationToEndOfSegment(p, theseAxons{i}, pos{i}{j} + [1 1 1], dir{i}{j}, 200, true);
            end
        end
    end
end
