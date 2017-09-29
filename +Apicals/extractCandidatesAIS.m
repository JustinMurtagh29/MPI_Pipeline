% extractCandidateAIS
% script to filter with directionality and spine head density the
% state without the outgrown soma to collect all AIS candidates
% -------------------------------------------------------------------------
% Author: Matej Zecevic <matej.zecevic@brain.mpg.de>


%% load states

% load outgrown soma removed state
% --------------------
load('/home/zecevicm/Desktop/connectomics_git/L4_apicalDendrites/getClassificationOfComponents/20170913_agglosNewFormatAfterSomaDeletion.mat');
% convert to old format
agglosSomaCut = connectEM.transformAggloNewOldRepr(newAgglos);
clear newAgglos;
% get meta info (points)
metaPath = '/mnt/gaba/gaba/u/mberning/results/pipeline/20170217_ROI/segmentMeta.mat';
sM = load(metaPath);
pointsSomaCut = cellfun(@(x) sM.point(:,x)',agglosSomaCut,'uni',0);
% load parameter
load('/mnt/gaba/gaba/u/mberning/results/pipeline/20170217_ROI/allParameterWithSynapses.mat');

% apply directionality
% ----------------------------------------
% load meta information from apicals script on original file
% agglos (all original agglos), rInd_x (indices after each filter),
% meanFPC (mean first pc, directionality)
load('/home/zecevicm/Desktop/connectomics_git/L4_apicalDendrites/getClassificationOfComponents/20170919_apicalsFromDend02_Meta.mat');

%% filter using directionality (var explained & mean FPC deviation)

tic;
% calculate PCA and var explained for each agglo
pcaMain = zeros(3, length(pointsSomaCut)); 
varExplained = zeros(length(pointsSomaCut),3);
for a=1:length(pointsSomaCut)
    D = pointsSomaCut{a};
    [coeff,score,latent,tsquared,explained,mu] = pca(D);
    pcaMain(:,a) = coeff(:,1);
    varExplained(a,:) = explained(1);
end
varExpFirst = varExplained(:,1);
disp('finished applying PCA to data.');

% kick out all agglos that are below threshold (don't explain enough
% variance) and get mean first principal component
varThr = 0.95;
varThr = varThr * 100;
pcaMain = pcaMain(:,transpose(varExpFirst >= varThr));
filteredAgglos = agglosSomaCut(varExpFirst >= varThr);
filteredPoints = pointsSomaCut(varExpFirst >= varThr);
meanFPCSomaCut = mean(pcaMain, 2);
% create new list of indices for all of somacut and then filter 
rIndSomaCut = 1:numel(agglosSomaCut);
rIndSC_1 = rIndSomaCut(varExpFirst >= varThr);

% filter that don't align with main axis
% use scalar product (angle)
spThr = 0.95;
sps = transpose(transpose(meanFPC) * pcaMain);
filteredAgglos = filteredAgglos(sps >= spThr);
filteredPoints = filteredPoints(sps >= spThr);
rIndSC_2 = rIndSC_1(sps >= spThr);

f1 = length(filteredAgglos) / length(agglos);
fprintf('filtered agglos, fraction of total: %f (%d/%d)\n', f1, length(filteredAgglos), length(agglos));

disp('finished execution of filterMainAxis.');
toc; disp('---------------------------------------');

% plot: histogram for distribution of distances (meanFPC from observed)
h = histogram(sps);
title([ num2str(length(sps)) ' agglos out of ' num2str(length(agglos)) ' - Bin width: ' num2str(h.BinWidth)]);
xlabel('distance from mean FPC'); ylabel('count');

%% take out all irrelevant straight things 
% smaller agglos etc. (all with <= e.g. 10 segments out)

l = cellfun(@length, filteredAgglos, 'un', 0);
rIndSC_3 = rIndSC_2(vertcat(l{:}) > 10);

fprintf('>10 segments, fraction of total: %d/%d \n', length(rIndSC_3), length(agglosSomaCut));
disp('---------------------------------------');

%% filter all that touch one border at least

% tolerance distance from border
tol = 300; %250; %200;

% only x boundaries (pia direction)
xbounds = p.bbox(1,:);

% filter out agglos that touch the border in x direction
mask = zeros(length(agglosSomaCut(rIndSC_3)), 1);
points = pointsSomaCut(rIndSC_3);

for i=1:length(agglosSomaCut(rIndSC_3))
    for j=1:size(points{i},1)
        point = points{i}(j);   % only consider x
        % if border is touched at one of the places within an agglo
        if point <= (xbounds(1) + tol)
            mask(i) = 1;
        elseif point >= (xbounds(2) - tol)
            mask(i) = 1;
        end
    end
end

rIndSC_4 = rIndSC_3(logical(mask));

fprintf('at least one x contact with %d tolerance, fraction of total: %d/%d \n', tol, length(rIndSC_4), length(agglosSomaCut));
disp('---------------------------------------');

%% get spine head counts and lengths

% spine head attachment
spineHeads = load('/home/zecevicm/Desktop/connectomics_git/L4_apicalDendrites/getClassificationOfComponents/20170919_spineheahOutgrownSomaCut.mat');
% path lengths in nm
dendLens = load('/home/zecevicm/Desktop/connectomics_git/L4_apicalDendrites/getClassificationOfComponents/20170919_pathlengthsOutgrownSomaCut.mat');
dendLens = dendLens.lens / 1E3; % in um

% get spine head counts
[ shCounts ] = Apicals.collectSHcounts( agglosSomaCut, spineHeads );

% collect spine head densities for the filtered Agglos
shDensity = shCounts ./ dendLens;
% replace NaN with zero spine density
shDensity(isnan(shDensity)) = 0; % zero spines

% plot: histogram for spine head densities
histogram(shDensity(rIndSC_4), 10);
xlabel('Spine density: # spines per um'); ylabel('# dendrite agglos');
title(['Spine Head density for ' num2str(length(rIndSC_4)) ' agglos filtered']);

%% filter with spine head density to get AIS candidates

shdThreshold = 0.15;
rIndSC_5 = rIndSC_4(shDensity(rIndSC_4) < shdThreshold);

% plot: histogram for spine head densities
histogram(shDensity(rIndSC_5), 10);
xlabel('Spine density: # spines per um'); ylabel('# dendrite agglos');
title(['Spine Head density for ' num2str(length(rIndSC_5)) ' agglos filtered']);

%% write

indices = rIndSC_5;

agglosToWrite =  agglosSomaCut(indices);
writePath = '/home/zecevicm/Desktop/connectomics_git/L4_apicalDendrites/20170922_candidatesAIS_174withSHD<015_withMeta.nml';
distthr = 3000;

tic;
% collect some meta to write out about the agglos
strings = arrayfun(@(x) strcat('skel_',num2str(x,'%.2i - ')), 1:numel(agglosToWrite),'uni',0);

% density, length, counts
aTW_shCounts = shCounts(indices); aTW_lens = dendLens(indices);
aTW_shDensity = shDensity(indices);
for i=1:length(agglosToWrite)
    strings(i) = strcat(strings(i), num2str(aTW_shDensity(i), 'shd: %2.4f - '), ...
        num2str(aTW_lens(i), 'length: %2.f - '), num2str(aTW_shCounts(i), 'shcount: %2.f'));
end

% actually write
connectEM.generateSkeletonFromNodes(writePath,...
    cellfun(@(x) sM.point(:,x)', agglosToWrite,'uni',0), ...
    strings,[],[],[],distthr); 
fprintf('created %d skeletons from filtered agglos.\n', length(agglosToWrite)); toc;
