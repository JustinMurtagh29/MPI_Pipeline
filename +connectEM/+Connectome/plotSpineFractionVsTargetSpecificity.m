% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';

connFile = fullfile( ...
    '/home/amotta/Desktop', ...
    'connectome_axons_18_a_with_den_meta.mat');

minSynPre = 10;
threshCount = 5;

info = Util.runInfo();

%% loading data
conn = load(connFile);

%% preparing data
axonMask = (conn.axonMeta.synCount >= minSynPre);

% target specificities
specificities = conn.classConnectome ./ sum(conn.classConnectome, 2);
specificities = specificities(axonMask, :);

% spine synapse fractions
spineFracs = conn.axonMeta.spineSynCount ./ conn.axonMeta.synCount;
spineFracs = spineFracs(axonMask, :);

classNames = conn.denClasses;
classCount = numel(classNames);

classNames{strcmpi(classNames, 'whole cells')} = 'WCs';
classNames{strcmpi(classNames, 'apical dendrites')} = 'ADs';
classNames{strcmpi(classNames, 'smooth dendrite')} = 'SDs';

%% splot spine fraction histogram as function of specificity threshold
fig = figure();

for curIdxC = 1:classCount
    curName = classNames{curIdxC};
    
    curSpecs = specificities(:, curIdxC);
    curSpecThreshs = linspace(0, max(curSpecs), threshCount + 1);
    curSpecThreshs(end) = [];
    
    for curIdxT = 1:threshCount
        curSpecThresh = curSpecThreshs(curIdxT);
        
        % filter axons and get spine fractions
        curMask = (curSpecs >= curSpecThresh);
        curSpineFracs = spineFracs(curMask);
        
        
        % plot histogram
        curAx = subplot( ...
            threshCount, classCount, ...
            curIdxC + (curIdxT - 1) * classCount);
        histogram(curAx, curSpineFracs, linspace(0, 1, 21));
        
        hold(curAx, 'on');
        
        % plot median
        curSpineFracMedian = median(curSpineFracs);
        plot( ...
            curAx, ...
           [curSpineFracMedian, curSpineFracMedian], ...
            curAx.YAxis.Limits, 'k--');
        
        xlim(curAx, [0, 1]);
        xticks(curAx, []);
        xticklabels(curAx, {});
        
        ylabel({ ...
            'Axons with'; sprintf( ...
            'S(%s) ≥ %.2f', curName, curSpecThresh)});
    end
    
    if curIdxT == threshCount && curIdxC == 1
        % panel at the bottom left
        xticks(curAx, [0, 1]);
        xticklabels(curAx, {'0', '1'});
        xlabel('Spine synapse fraction');
    end
end

annotation( ...
    fig, ...
    'textbox', [0, 0.9, 1, 0.1], ...
	'String', { ...
        'Spine fraction vs. target specificity'; ...
        info.git_repos{1}.hash}, ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center');

fig.Position(3:4) = [1640, 1090];
