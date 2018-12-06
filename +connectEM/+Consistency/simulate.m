% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
numSyn = 2;
numSteps = 2;
numRuns = 1000;

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
for curRun = 1:numRuns
    curHasStdp = rand(1) < stdpProb;
    wMat(:, 1, curRun) = min(1, lognrnd(mu, sigma, [2, 1]));

    for curStep = 2:numSteps
        curOldW = wMat(:, curStep - 1, curRun);
        curDeltaT = 2 * tau * (2 * rand(1) - 1);

        % HACK(amotta): Correlate strength and 
        % curDeltaT = curDeltaT + 2 * sum(curOldW);

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
figure;

subplot(2, 1, 1);
hold on;
plot(wMat(1, :, 1));
plot(wMat(2, :, 1));
ylim([0, 1]);

subplot(2, 1, 2);
plot(std(wMat(:, :, 1), 0, 1) ./ mean(wMat(:, :, 1), 1));
ylim([0, sqrt(2)]);
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

figure;
hold on;

histogram(cvs(:, 1, :), ...
    'BinEdges', curEdges, ...
    'Normalization', 'probability', ...
    'DisplayStyle', 'stairs', 'LineWidth', 2);
histogram(cvs(:, end, :), ...
    'BinEdges', curEdges, ...
    'Normalization', 'probability', ...
    'DisplayStyle', 'stairs', 'LineWidth', 2);
histogram(cvs(:, end, ltpMask), ...
    'BinEdges', curEdges, ...
    'Normalization', 'probability', ...
    'DisplayStyle', 'stairs', 'LineWidth', 2);
histogram(cvs(:, end, ltdMask), ...
    'BinEdges', curEdges, ...
    'Normalization', 'probability', ...
    'DisplayStyle', 'stairs', 'LineWidth', 2);

%%
asis = 10 .^ mean(log10(wMat), 1);
curEdges = linspace(0, 1, 21);

figure;
hold on;

histogram(asis(:, 1, :), ...
    'BinEdges', curEdges, ...
    'DisplayStyle', 'stairs', 'LineWidth', 2);
histogram(asis(:, end, :), ...
    'BinEdges', curEdges, ...
    'DisplayStyle', 'stairs', 'LineWidth', 2);

%{
curBefore = histcounts(asis(:, 1, :), curEdges);
curLtp = histcounts(asis(:, 1, ltpMask), curEdges);
curLtpEst = 0.125 * curBefore;

histogram( ...
    'BinEdges', curEdges, 'BinCounts', max(curBefore - curLtp, 0), ...
    'DisplayStyle', 'stairs', 'LineWidth', 2)

histogram(asis(:, 1, ltpMask), ...
    'BinEdges', curEdges, ...
    'DisplayStyle', 'stairs', 'LineWidth', 2);
histogram(asis(:, end, ltdMask), ...
    'BinEdges', curEdges, ...
    'DisplayStyle', 'stairs', 'LineWidth', 2);

curKnots = unique(-asis(:, [1, end], :));
curA = histcounts(-asis(:, 1, :), curKnots, 'Normalization', 'cdf');
curB = histcounts(-asis(:, end, :), curKnots, 'Normalization', 'cdf');

[~, curThresh] = max(curB - curA);
curThresh = -curKnots(curThresh);

curLtpFrac = ...
    mean(asis(:, end, :) > curThresh) ...
  - mean(asis(:, 1, :) > curThresh);
curLtpFrac = 0.125;

curBefore = histcounts(asis(:, 1, :), curEdges);
curAfter = histcounts(asis(:, end, :), curEdges);
curLtp = histcounts(asis(:, 1, ltpMask), curEdges);
% curLtd = curAfter - (curBefore - curLtpFrac * curBefore);
curLtd = curAfter - (curBefore - curLtp);

histogram('BinEdges', curEdges, 'BinCounts', max(curLtd, 0))
%}
    
%%
%{
curEdges = linspace(-15, 0, 31);
curPBefore = sum(log(logncdf(wMat(:, 1, :), mu, sigma)), 1);
curPAfter = sum(log(logncdf(wMat(:, end, :), mu, sigma)), 1);

figure; hold on;
histogram(curPBefore, 'BinEdges', curEdges, 'DisplayStyle', 'stairs', 'LineWidth', 2);
histogram(curPAfter, 'BinEdges', curEdges, 'DisplayStyle', 'stairs', 'LineWidth', 2);

curAx = gca;
curAx.YScale = 'log';
%}

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

%{
annotation( ...
    curFig, ...
    'textbox', [0, 0.9, 1, 0.1], ...
    'String', { ...
    info.filename; info.git_repos{1}.hash; ...
    curConfig.title; curCtrlConfig.title}, ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center');
%}


%%

curEdges = linspace(0, 1, 21);
eW = wMat(:, end, :);
wCv = std(eW, 0, 1) ./ mean(eW, 1);
eMask = wCv < 0.4;

figure;
hold on;
histogram(eW, 'BinEdges', curEdges, 'Normalization', 'probability');
histogram(eW(:, :, eMask), 'BinEdges', curEdges, 'Normalization', 'probability');

%% Simulate
%{
curInitWeights = wMat(:, end, :);

curInitCv = std(curInitWeights, 0, 1) ./ mean(curInitWeights, 1);
curInitWeights = curInitWeights(:, :, curInitCv < 0.4);

rng(0);
for curRun = 1:numRuns
    curRandIds = randperm(numel(curInitWeights), 2);
    wMat(:, 1, curRun) = curInitWeights(curRandIds);

    for curStep = 2:numSteps
        curOldW = wMat(:, curStep - 1, curRun);
        curDeltaT = 2 * tau * (2 * rand(1) - 1);

        % HACK(amotta): Correlate strength and 
        % curDeltaT = curDeltaT + 2 * sum(curOldW);

        if curDeltaT > 0
            curDeltaW = +nuPos * aPos(curOldW);
        else
            curDeltaW = -nuNeg * aNeg(curOldW);
        end

        % Make synapses unreliable
        curErr = rand(numSyn, 1) < errProb;
        curDeltaW = (1 - curErr) .* curDeltaW;

        wMat(:, curStep, curRun) = curOldW + curDeltaW;
    end
end
%}
