% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
numSyn = 2;
numSteps = 2;
numRuns = 10000;

tau = 10; % ms

wMax = 1;
nuPos = 0.5;
aPos = @(w) (wMax - w) * nuPos;

nuNeg = 0.5;
aNeg = @(w) w * nuNeg;

wMat = nan(numSyn, numSteps, numRuns);

errProb = 0;
stdpProb = 0.25;

mu = -2;
sigma = 1;

%% Simulate
rng(0);

% Normal initialization
wMat(:, 1, :) = min(1, lognrnd(mu, sigma, [2, 1, numRuns]));

for curRun = 1:numRuns
    curHasStdp = rand(1) < stdpProb;

    for curStep = 2:numSteps
        curOldW = wMat(:, curStep - 1, curRun);
        curDeltaT = 2 * tau * (2 * rand(1) - 1);

        % HACK(amotta): Correlate strength and 
        % curDeltaT = curDeltaT + 4 * sum(curOldW);

        if curDeltaT > 0
            curDeltaW = +aPos(curOldW) * exp(-curDeltaT / tau);
        else
            curDeltaW = -aNeg(curOldW) * exp(+curDeltaT / tau);
        end

        % Make synapses unreliable
        curErr = rand(numSyn, 1) < errProb;
        curDeltaW = (1 - curErr) .* curDeltaW;
        
        curDeltaW = curDeltaW * curHasStdp;

        wMat(:, curStep, curRun) = curOldW + curDeltaW;
    end
end

ltpMask = all(wMat(:, end, :) > wMat(:, 1, :), 1);
ltdMask = all(wMat(:, end, :) < wMat(:, 1, :), 1);

%%
%{
curFig = figure;
curFig.Color = 'white';
curFig.Position(3:4) = [355, 319];

subplot(2, 1, 1);
hold on;
plot(wMat(1, :, 1), 'LineWidth', 2);
plot(wMat(2, :, 1), 'LineWidth', 2);

ylim([0, 1]);
ylabel('Synaptic weights');
set(gca, 'TickDir', 'out');

subplot(2, 1, 2);
hold on;
plot(std(wMat(:, :, 1), 0, 1) ./ mean(wMat(:, :, 1), 1), 'black', 'LineWidth', 2);
ylim([0, sqrt(2)]);
ylabel('CV');
set(gca, 'TickDir', 'out');
%}

%%
%{
cvs = std(wMat, 0, 1) ./ mean(wMat, 1);
cvMeans = reshape(mean(cvs, 3), [], 1);
cvStd = reshape(std(cvs, 0, 3), [], 1);

figure;
hold on;
plot(cvMeans, 'black');
plot(cvMeans + cvStd, 'k--');
plot(cvMeans - cvStd, 'k--');
ylim([0, sqrt(2)]);
%}

%%
cvs = std(wMat, 0, 1) ./ mean(wMat, 1);
curEdges = linspace(0, 1.5, 16);

curFig = figure();
curFig.Color = 'white';
curFig.Position(3:4) = [205, 165];

curAx = axes(curFig);
axis(curAx, 'square');
curAx.TickDir = 'out';
hold on;

histogram(cvs(:, 1, :), ...
    'BinEdges', curEdges, ...
    'DisplayStyle', 'stairs', 'LineWidth', 2);
histogram(cvs(:, end, :), ...
    'BinEdges', curEdges, ...
    'DisplayStyle', 'stairs', 'LineWidth', 2);
histogram(cvs(:, end, ltpMask), ...
    'BinEdges', curEdges, ...
	'DisplayStyle', 'stairs', 'LineWidth', 2);

xlim(curAx, [0, 1.5]);
xlabel(curAx, 'Coefficient of variation');
ylabel(curAx, 'Connections');

curKnots = unique(cvs(:, [1, end], :));
curBefore = histcounts(cvs(:, 1, :), curKnots, 'Normalization', 'cdf');
curAfter = histcounts(cvs(:, end, :), curKnots, 'Normalization', 'cdf');

[~, curThreshIdx] = max(curAfter - curBefore);
curThresh = curKnots(curThreshIdx);
curAfter(curThreshIdx) - curBefore(curThreshIdx) %#ok

%%
% Biased initialization
curInitWeights = wMat(:, end, :);
wMat(:, 1, :) = curInitWeights(randi(numel(curInitWeights), [2, 1, numRuns]));

%%
curSaSdT = table;
curSaSdT.areas = transpose(sort(squeeze(wMat(:, end, :)), 1));
curSaSdT.cv = std(curSaSdT.areas, 0, 2) ./ mean(curSaSdT.areas, 2);
curSaSdT.avgLogAreas = mean(log10(curSaSdT.areas), 2);

curCtrlT = table;
curCtrlT.areas = transpose(sort(squeeze(wMat(:, 1, :)), 1));
curCtrlT.cv = std(curCtrlT.areas, 0, 2) ./ mean(curCtrlT.areas, 2);
curCtrlT.avgLogAreas = mean(log10(curCtrlT.areas), 2);

% Density difference map
curLimX = [-eps, 1.5];
curLimY = [-1.5, 0.5];

curTicksX = linspace(curLimX(1), curLimX(2), 4);
curTicksY = linspace(curLimY(1), curLimY(2), 5);

curImSize = [301, 301];

[curImGridY, curImGridX] = ndgrid( ...
    linspace(curLimY(1), curLimY(2), curImSize(1)), ...
    linspace(curLimX(1), curLimX(2), curImSize(2)));
curImGrid = cat(2, curImGridY(:), curImGridX(:));

curSaSdImg = horzcat( ...
    curSaSdT.avgLogAreas, curSaSdT.cv);
curMask = ...
    curLimY(1) <= curSaSdImg(:, 1) ...
    & curSaSdImg(:, 1) <= curLimY(2) ...
    & curLimX(1) <= curSaSdImg(:, 2) ...
    & curSaSdImg(:, 2) <= curLimX(2);

curSaSdImg = ksdensity( ...
    curSaSdImg(curMask, :), curImGrid, ...
    'Support', transpose(cat(1, curLimY, curLimX)), ...
    'BoundaryCorrection', 'reflection');
curSaSdImg = mean(curMask) * curSaSdImg / sum(curSaSdImg(:));
curSaSdImg = reshape(curSaSdImg, curImSize);

curCtrlImg = horzcat( ...
    curCtrlT.avgLogAreas, curCtrlT.cv);
curMask = ...
    curLimY(1) <= curCtrlImg(:, 1) ...
    & curCtrlImg(:, 1) <= curLimY(2) ...
    & curLimX(1) <= curCtrlImg(:, 2) ...
    & curCtrlImg(:, 2) <= curLimX(2);

curCtrlImg = ksdensity( ...
    curCtrlImg(curMask, :), curImGrid, ...
    'Support', transpose(cat(1, curLimY, curLimX)), ...
    'BoundaryCorrection', 'reflection');
curCtrlImg = mean(curMask) * curCtrlImg / sum(curCtrlImg(:));
curCtrlImg = reshape(curCtrlImg, curImSize);

curMax = max(max(curSaSdImg(:)), max(curCtrlImg(:)));
curDiffImg = curSaSdImg - curCtrlImg;
curMaxDiff = max(abs(curDiffImg(:)));

curFig = figure();
curFig.Color = 'white';
curFig.Position(3:4) = [360, 1020];

curAx = subplot(3, 1, 1);
image(curAx, ...
    uint8(double(intmax('uint8')) ...
    * curSaSdImg / curMax));
colormap(curAx, jet(256));

curBar = colorbar('peer', curAx);
curBar.Ticks = curBar.Limits;
curBar.TickLabels = {'0', sprintf('%.3g', curMax)};

curAx = subplot(3, 1, 2);
image(curAx, ...
    uint8(double(intmax('uint8')) ...
    * curCtrlImg / curMax));
colormap(curAx, jet(256));

curBar = colorbar('peer', curAx);
curBar.Ticks = curBar.Limits;
curBar.TickLabels = {'0', sprintf('%.3g', curMax)};

curAx = subplot(3, 1, 3);
image(curAx, ...
    uint8(double(intmax('uint8')) ...
    * (1 + curDiffImg / curMaxDiff) / 2));
colormap(curAx, jet(256));
curBar = colorbar('peer', curAx);
curBar.Label.String = { ...
    'Fraction of observe pairs'; ...
    'relative to null model'};
curBar.Ticks = [ ...
    curBar.Limits(1), ...
    mean(curBar.Limits), ...
    curBar.Limits(end)];
curBar.TickLabels = { ...
    sprintf('%.3g', -curMaxDiff), '0', ...
    sprintf('%.3g', +curMaxDiff)};

set( ...
    findobj(curFig.Children, 'Type', 'ColorBar'), ...
    'Location', 'EastOutside', ...
    'TickDirection', 'out', ...
    'Box', 'off');

[~, curTickIdsX] = ismember(curTicksX, ...
    linspace(curLimX(1), curLimX(2), curImSize(2)));
curTickLabelsX = arrayfun( ...
    @num2str, curTicksX, 'UniformOutput', false);
[~, curTickIdsY] = ismember(curTicksY, ...
    linspace(curLimY(1), curLimY(2), curImSize(1)));
curTickLabelsY = arrayfun( ...
    @num2str, curTicksY, 'UniformOutput', false);

curAxes = reshape(flip(findobj(curFig, 'Type', 'Axes')), 1, []);
arrayfun(@(ax) hold(ax, 'on'), curAxes);

set(curAxes, ...
    'Box', 'off', ...
    'TickDir', 'out', ...
    'YDir', 'normal', ...
    'YTick', curTickIdsY, ...
    'YTickLabels', curTickLabelsY, ...
    'XTick', [], ...
    'PlotBoxAspectRatio', [1, 1, 1], ...
    'DataAspectRatioMode', 'auto');
set(curAxes(end), ...
    'XTick', curTickIdsX, ...
    'XTickLabels', curTickLabelsX);

arrayfun(@(ax) ylabel(ax, ...
    'Average log_{10}(ASI area [µm²])'), curAxes);
xlabel(curAxes(end), 'Coefficient of variation');

%% LTP
numSteps = 5;

ltpMat = [0; 0.75];
for curStep = 1:(numSteps - 1)
    ltpMat(:, curStep + 1) = ltpMat(:, curStep) + 0.5 * (1 - ltpMat(:, curStep));
end

tVec = repelem(0:numSteps, 2);
tVec = tVec(2:(end - 1));
ltpMat = repelem(ltpMat, 1, 2);

curFig = figure;
curFig.Color = 'white';
curFig.Position(3:4) = [250, 90];
hold on;
plot(tVec, ltpMat(1, :), 'black', 'LineWidth', 2);
plot(tVec, ltpMat(2, :), 'black', 'LineWidth', 2);

set(gca, 'TickDir', 'out');
xlim(tVec([1, end]));

%% LTP
numSteps = 5;

ltdMat = [1; 0.25];
for curStep = 1:(numSteps - 1)
    ltdMat(:, curStep + 1) = ltdMat(:, curStep) - 0.5 * ltdMat(:, curStep);
end

tVec = repelem(0:numSteps, 2);
tVec = tVec(2:(end - 1));
ltdMat = repelem(ltdMat, 1, 2);

curFig = figure;
curFig.Color = 'white';
curFig.Position(3:4) = [250, 90];
hold on;
plot(tVec, ltdMat(1, :), 'black', 'LineWidth', 2);
plot(tVec, ltdMat(2, :), 'black', 'LineWidth', 2);

set(gca, 'TickDir', 'out');
xlim(tVec([1, end]));

