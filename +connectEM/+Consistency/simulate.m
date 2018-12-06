% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
numSyn = 2;
numSteps = 400;
numRuns = 1000;

tau = 10; % ms

wMax = 1;
nuPos = 0.2;
aPos = @(w) (wMax - w) * nuPos;

nuNeg = 0.2;
aNeg = @(w) w * nuNeg;

wMat = nan(numSyn, numSteps, numRuns);

errProb = 0.25;

%% Simulate
rng(0);
for curRun = 1:numRuns
    wMat(:, 1, curRun) = min(1, lognrnd(-1.65, 1, [2, 1]));

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

%%
figure;

subplot(2, 1, 1);
hold on;
plot(wMat(1, :, 1));
plot(wMat(2, :, 1));
ylim([0, 1]);

subplot(2, 1, 2);
plot(std(wMat(:, :, 1), 0, 1) ./ mean(wMat(:, :, 1), 1));
ylim([0, sqrt(2)]);

%%
cvs = std(wMat, 0, 1) ./ mean(wMat, 1);
cvMeans = reshape(mean(cvs, 3), [], 1);
cvStd = reshape(std(cvs, 0, 3), [], 1);

figure;
hold on;
plot(cvMeans, 'black');
plot(cvMeans + cvStd, 'k--');
plot(cvMeans - cvStd, 'k--');
ylim([0, sqrt(2)]);

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

%%
asis = 10 .^ mean(log10(wMat), 1);
curEdges = linspace(0, 1, 21);

figure;
hold on;

histogram(asis(:, 1, :), ...
    'BinEdges', curEdges, ...
    'Normalization', 'probability', ...
    'DisplayStyle', 'stairs', 'LineWidth', 2);
histogram(asis(:, end, :), ...
    'BinEdges', curEdges, ...
    'Normalization', 'probability', ...
    'DisplayStyle', 'stairs', 'LineWidth', 2);

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
