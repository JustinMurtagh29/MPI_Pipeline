% Note
%   All data in variables prefixed by `mish` are from
%   Mishchenko et al. (2010) Neuron
%
% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
info = Util.runInfo();
Util.showRunInfo(info);

%% Data from panel 3d of Mishchenko et al. (2010) Neuron
% Extracted from PDF using Inkskape...
mishData = [ ...
    395.998, 651.110;
    347.881, 636.625;
    358.979, 627.814;
    409.723, 633.151;
    466.649, 639.250;
    434.967, 621.714;
    401.335, 616.717;
    375.751, 609.517;
    383.801, 607.145;
    366.603, 607.567;
    369.653, 607.145;
    373.126, 606.382;
    378.803, 606.043;
    398.709, 607.907;
    355.929, 605.619;
    358.217, 604.518;
    365.417, 603.757;
    373.126, 604.518;
    373.466, 601.807;
    376.938, 599.182;
    390.661, 601.045;
    423.107, 599.519;
    393.373, 596.894;
    434.205, 594.606;
    366.603, 596.894;
    354.405, 586.219;
    370.078, 592.319;
    395.659, 590.033;
    375.414, 582.323;
    392.949, 582.323;
    334.158, 583.086;
    360.505, 570.887];

mishDiam = 4.577;
mishLims = [0, 7];
mishRangeX = 318.826 + [0, 152.060];
mishRangeY = 537.255 + [0, 152.061];

mishData = mishData + mishDiam / 2;
mishData = mishData - [mishRangeX(1), mishRangeY(1)];
mishData = mishData ./ [diff(mishRangeX), diff(mishRangeY)];
mishData = mishLims(1) + diff(mishLims) .* mishData;

% Dendrite length. Mean their two measurements. Assumed to be constant
% because it primarily reflects extents of reconstructed volumes.
mishDendLen = (3.45 + 3.64) / 2;

mishRsq = 0.12;

%% Reproducing panel 3d of Mishchenko et al. (2010) Neuron
clear cur*;

curFit = fitlm(mishData(:, 2), mishData(:, 1));

curFig = figure();
curAx = axes(curFig);
hold(curAx, 'on');

plot(curAx, mishLims, mishLims, 'k--');
scatter(curAx, mishData(:, 1), mishData(:, 2), 'ko');
curFitPlot = plot(curAx, curFit.predict(mishLims(:)), mishLims(:), 'k');

xlim(curAx, mishLims); ylim(curAx, mishLims);
xlabel(curAx, 'Actual density of synapses (µm^{-1})');
ylabel(curAx, 'Predicted density of synapses (µm^{-1})');
axis(curAx, 'square');

curLeg = sprintf('Linear regression (r² = %g)', curFit.Rsquared.Ordinary);
curLeg = legend(curFitPlot, curLeg, 'Location', 'NorthWest'); %#ok

connectEM.Figure.config(curFig, info);
curAx.Title.String{end + 1} = 'Mishchenko et al. (2010) Neuron, Panel 3d';

%% Forward model
% The idea of the forward model is the following: Let's assume that the
% second equation, and the data shown in panel 3d of Mishchenko et al.
% (2010) Neuron correctly predict the dendritic synapse density.
%   If we now assume that dendrites establish synapses following a Poisson
% distribution with the predicted synapse density, we can estimate the
% amount of variance that is due to the dendrite stretches being limited to
% around 3.5 µm in length.
%   This also yields a p-value for the reported r² value.
clear cur*;
rng(0);

% NOTE(amotta): Build synapse counts by sampling according to axonal
% density. A Poisson distribution is used.
mishPredSynDensity = mishData(:, 2);
synDensity = poissrnd(mishDendLen .* mishPredSynDensity) ./ mishDendLen;

curFit = fitlm(mishData(:, 2), synDensity);

curFig = figure();
curAx = axes(curFig);
hold(curAx, 'on');

plot(curAx, mishLims, mishLims, 'k--');
scatter(curAx, synDensity, mishData(:, 2), 'ko');
curFitPlot = plot(curAx, curFit.predict(mishLims(:)), mishLims(:), 'k');

% xlim(curAx, mishLims); ylim(curAx, mishLims);
xlabel(curAx, 'Synthetic density of synapses (µm^{-1})');
ylabel(curAx, 'Predicted density of synapses (µm^{-1})');
axis(curAx, 'square');

curLeg = sprintf('Linear regression (r² = %g)', curFit.Rsquared.Ordinary);
curLeg = legend(curFitPlot, curLeg, 'Location', 'NorthWest');

connectEM.Figure.config(curFig, info);
curAx.Title.String{end + 1} = 'Forward model';

%% Estimate p-value for reported r-squared value
clear cur*;
rng(0);

curN = 10000;
curRsq = nan(curN, 1);
for curIdx = 1:curN
    curSynDens = mishDendLen .* mishPredSynDensity;
    curSynDens = poissrnd(curSynDens) ./ mishDendLen;
    
    curFit = fitlm(mishData(:, 2), curSynDens);
    curRsq(curIdx) = curFit.Rsquared.Ordinary;
end

curPVal = mean(curRsq <= mishRsq) %#ok
curMedianRsq = median(curRsq) %#ok
curMeanRsq = mean(curRsq) %#ok
