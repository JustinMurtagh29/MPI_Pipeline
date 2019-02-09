% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
nmlDir = fullfile(fileparts(mfilename('fullpath')), 'annotations');

% HACKHACKHACKHACK(amotta): Fix this!
asiFile = '/mnt/mpibr/data/Personal/mottaa/L4/2019-02-08-ASI-Area-Calibration/asiT-with-areas.mat';

methodNames = { ...
    'Sum of Border Areas', ...
    'Helmstaedter et al. 2013 Nature', ...
    'de Vivo et al. 2017 Science', ...
    'Isosurface Area', ...
    'Alpha Shape Area'};

info = Util.runInfo();
Util.showRunInfo(info);

%% Load data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

asiT = load(asiFile);
asiT = asiT.asiT;

%% Load and prepare calibration data
calibT = connectEM.Consistency.Calibration.loadAsiAreaNmls(param, nmlDir);

[~, autoAreas] = ismember(calibT.asiId, asiT.id);
autoAreas = [asiT.area(autoAreas), asiT.areas(autoAreas, :)];
assert(size(autoAreas, 2) == numel(methodNames));

%% Show calibration data
for curId = 1:numel(methodNames)
    curX = log(autoAreas(:, curId));
    curY = log(calibT.area);
    
    curCoeff = ones(size(curX)) \ (curY - curX);
    curRmse = sqrt(mean((curY - (curX + curCoeff)) .^ 2));
    
    curFig = figure();
    curAx = axes(curFig);
    hold(curAx, 'on');
    
    scatter(curAx, curX, curY, 80, '.');
    curLim = prctile(cat(2, xlim(curAx), ylim(curAx)), [0, 100]);
    plot(curLim, curLim + curCoeff);
    
    plot(curAx, ...
        curLim, curLim, ...
        'Color', 'black', 'LineStyle', '--');
    
    curFig.Color = 'white';
    curFig.Position(3:4) = 330;
    
    axis(curAx, 'square');
    set(curAx, 'XLim', curLim, 'YLim', curLim, 'TickDir', 'out');
    
    xlabel(curAx, 'Automated log_{e}(ASI area [µm²])');
    ylabel(curAx, 'Tracing-based log_{e}(ASI area [µm²])');
    
    curLeg = { ...
        sprintf('%d calibration points', height(calibT)), ...
        sprintf( ...
            'Factor %.3g (%.3g RMSE)', exp(curCoeff), curRmse)};
    curLeg = legend(curAx, curLeg, 'Location', 'NorthWest');
    curLeg.Box = 'off';
    
    title(curAx, ...
        {info.filename; info.git_repos{1}.hash; methodNames{curId}}, ...
        'FontWeight', 'normal', 'FontSize', 10);
end
