% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
numSyn = 2;
numSteps = 2;
numRuns = 10000;

tau = 10; % ms

wMax = 1;
nuPos = 0.25;
aPos = @(w) (wMax - w) * nuPos;

nuNeg = 0.25;
aNeg = @(w) w * nuNeg;

wMat = nan(numSyn, numSteps, numRuns);

errProb = 0;
stdpProb = 0.25;

%% Simulate
rng(0);
for curRun = 1:numRuns
    curHasStdp = rand(1) < stdpProb;
    wMat(:, 1, curRun) = min(1, lognrnd(-2, 1, [2, 1]));

    for curStep = 2:numSteps
        curOldW = wMat(:, curStep - 1, curRun);
        curDeltaT = 2 * tau * (2 * rand(1) - 1);

        % HACK(amotta): Correlate strength and 
        % curDeltaT = curDeltaT + 2 * sum(curOldW);

        if curDeltaT > 0
            curDeltaW = +aPos(curOldW);
        else
            curDeltaW = -aNeg(curOldW);
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

%% Simulate
%{
curInitWeights = wMat(:, end, :);

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
