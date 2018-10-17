load('/gaba/u/mberning/results/pipeline/20170217_ROI/allParameterWithSynapses.mat');
disp('Parameters loaded');
outputFolder = fullfile(p.saveFolder, 'aggloState');

%% load newest whole cell state
load(fullfile(outputFolder,'wholeCells_07.mat'))
[wcLUT,wcSegIds] = Superagglos.buildLUT(wholeCells);


%% load dend/axon state before any interaction
load(fullfile(outputFolder,'dendrites_03_v2.mat'),'dendrites','indBigDends');
load(fullfile(outputFolder,'axons_03.mat'),'axons','indBigAxons');
if islogical(indBigDends)
    indBigDends = find(indBigDends);
end
if islogical(indBigAxons)
    indBigAxons = find(indBigAxons);
end
[axonLUT,axonSegIds] = Superagglos.buildLUT(axons);
[dendriteLUT,dendriteSegIds] = Superagglos.buildLUT(dendrites);

%% make luts the same size
maxSegId = max([max(axonSegIds),max(dendriteSegIds),max(wcSegIds)]);
if numel(axonLUT) < maxSegId
    axonLUT(maxSegId) = 0;
end
if numel(dendriteLUT) < maxSegId
    dendriteLUT(maxSegId) = 0;
end
if numel(wcLUT) < maxSegId
    wcLUT(maxSegId) = 0;
end

%% calculate number of mergers and splits
mergers = zeros(numel(wholeCells),2);
splits = zeros(numel(wholeCells),2);
splitsBig = zeros(numel(wholeCells),2);

for w = 1:numel(wholeCells)
    disp(w)
    % overlap dend state with whole cells -> numel dends per whole cell -1
    % are the splits
    overlappingAgglos = unique(nonzeros(dendriteLUT(wcLUT == w)));
    if numel(overlappingAgglos) > 0
        splits(w,1) = splits(w,1) + numel(overlappingAgglos) - 1;
        splitsBig(w,1) = splitsBig(w,1) + numel(intersect(overlappingAgglos,indBigDends)) - 1;
    end
    % avoid pollution of counting by small axonic stuff which was already
    % taken into account by the addition of dendrites
    axonLUT(ismember(dendriteLUT,overlappingAgglos)) = 0;
    
    % now calculate the number of mergers by deleting edges of all
    % overlapping agglos that are connecting the correct with the incorrect
    % stuff and count the number of connected components
    for o = 1:numel(overlappingAgglos)
        ismem = ismember(dendrites(overlappingAgglos(o)).nodes(:,4),wholeCells(w).nodes(:,4));
        if any(~ismem)
            edgesToDelete = (sum(ismember(dendrites(overlappingAgglos(o)).edges,find(ismem)),2)==1);
            mergers(w,1) = mergers(w,1) + numel(Superagglos.removeEdgesFromAgglo(dendrites(overlappingAgglos(o)),edgesToDelete))-1;
        end
    end
    
    % overlap axon state with whole cells -> numel axons per whole cell -1
    % are the splits
    overlappingAgglos = unique(nonzeros(axonLUT(wcLUT == w)));
    if numel(overlappingAgglos) > 0
        splits(w,2) = splits(w,2) + numel(overlappingAgglos) - 1;
        splitsBig(w,2) = splitsBig(w,2) + numel(intersect(overlappingAgglos,indBigAxons)) - 1;
    end
    % now calculate the number of mergers by deleting edges of all
    % overlapping agglos that are connecting the correct with the incorrect
    % stuff and count the number of connected components
    for o = 1:numel(overlappingAgglos)
        ismem = ismember(axons(overlappingAgglos(o)).nodes(:,4),wholeCells(w).nodes(:,4));
        if any(~ismem)
            edgesToDelete = (sum(ismember(axons(overlappingAgglos(o)).edges,find(ismem)),2)==1);
            mergers(w,2) = mergers(w,2) + numel(Superagglos.removeEdgesFromAgglo(axons(overlappingAgglos(o)),edgesToDelete))-1;
        end
    end
end

fprintf('Total number of dendrite splits: %d\n',sum(splits(:,1)))
fprintf('Total number of axon splits: %d\n',sum(splits(:,2)))
fprintf('Total number of big dendrite splits: %d\n',sum(splitsBig(:,1)))
fprintf('Total number of big axon splits: %d\n',sum(splitsBig(:,2)))
fprintf('Total number of dendrite merger: %d\n',sum(mergers(:,1)))
fprintf('Total number of axon merger: %d\n',sum(mergers(:,2)))
    