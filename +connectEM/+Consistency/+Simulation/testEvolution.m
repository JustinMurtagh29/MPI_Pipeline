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

colors = get(groot, 'defaultAxesColorOrder');
colors = colors(1:size(methods, 1), :);

info = Util.runInfo();
Util.showRunInfo(info);

%% Plot
rng(4);
areas = 10 .^ normrnd(ccMeanLog, ccStdLog, [3, 2, size(methods, 1)]);

fig = figure();
ax = axes(fig);
hold(ax, 'on');

for curMethodIdx = 1:size(methods, 1)
    curMethod = methods{curMethodIdx, 1};
    curColor = colors(curMethodIdx, :);
    
    for curAreaIdx = 1:size(areas, 1)
        curAreas = ...
            connectEM.Consistency.Simulation.evolution( ...
                areas(curAreaIdx, :, curMethodIdx), 'method', curMethod);

        curX = std(curAreas, 0, 2) ./ mean(curAreas, 2);
        curY = log10(mean(curAreas, 2));
        
        scatter( ...
            ax, curX(1), curY(1), 'o', ...
            'MarkerEdgeColor', 'none', ...
            'MarkerFaceColor', curColor);
        plot(ax, curX, curY, 'Color', curColor);
    end
end

xlim(ax, limX);
ylim(ax, limY);

lines = findobj(ax, 'type', 'line');
set(lines, 'LineWidth', 2);

xlabel(ax, 'Coefficient of variation');
ylabel(ax, 'log_{10}(Average ASI area [µm²])');
axis(ax, 'square');

leg = flip(lines(1:size(areas, 1):end));
leg = legend(leg, methods(:, 2), 'Location', 'SouthOutside');
leg.Box = 'off';

fig.Position(3:4) = [320, 420];
connectEM.Figure.config(fig, info);
