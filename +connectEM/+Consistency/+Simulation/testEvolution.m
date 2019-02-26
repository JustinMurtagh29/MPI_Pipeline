% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
% Rough parameters for CC synapses
ccMeanLog = -0.707029;
ccStdLog = 0.298972;

limX = [0, 1.5];
limY = [-1.5, 0.5];

methods = { ...
    'ltdSubtractive', 'LTD. Linearly.'; ...
    'ltdDivisive',    'LTD. Exponentially.'; ...
    'ltdSaturated',   'LTD. Exponentially towards non-zero.'; ...
    'ltpSaturated',   'LTP. Exponentially towards maximum.'};

info = Util.runInfo();
Util.showRunInfo(info);

%% Plot
rng(0);
areas = 10 .^ normrnd(ccMeanLog, ccStdLog, [20, 2]);

fig = figure();
fig.Position(3:4) = [1330, 340];

for curMethodIdx = 1:size(methods, 1)
    subplot(1, size(methods, 1), curMethodIdx);
    hold('on');
    
    for curAreaIdx = 1:size(areas, 1)
        curAreas = ...
            connectEM.Consistency.Simulation.evolution( ...
                areas(curAreaIdx, :), 'method', methods{curMethodIdx, 1});

        curX = std(curAreas, 0, 2) ./ mean(curAreas, 2);
        curY = log10(mean(curAreas, 2));
        plot(curX, curY);
    end
    
    title(methods{curMethodIdx, 2});

    xlim(limX);
    ylim(limY);
end

lines = findobj(fig, 'type', 'line');
set(lines, 'LineWidth', 2);

ax = fig.Children(end);
xlabel(ax, 'Coefficient of variation');
ylabel(ax, 'log_{10}(Average ASI area [µm²])');

set( ...
    fig.Children, ...
    'PlotBoxAspectRatio', [1, 1, 1], ...
    'DataAspectRatioMode', 'auto');

connectEM.Figure.config(fig, info);
