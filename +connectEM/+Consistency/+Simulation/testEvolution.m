% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
% Rough parameters for CC synapses
ccMeanLog = -0.707029;
ccStdLog = 0.298972;

limX = [0, 1.5];
limY = [-1.5, 0.5];

curReds = autumn(2 + 1);
curGreens = summer(4 + 1);

methods = { ...
    'ltdLinearZero', 'LTD. Linearly towards zero.', curGreens(1, :); ...
    'ltdLinearNonZero', 'LTD. Linearly towards non-zero.', curGreens(2, :);
    'ltdExponentialZero', 'LTD. Exponential decay to zero.', curGreens(3, :); ...
    'ltdExponentialNonZero', 'LTD. Exponential decay to non-zero.', curGreens(4, :); ...
    'ltpExponentialFinite', 'LTP. Exponential grows to maximum.', curReds(1, :); ...
    'ltpLinearInf', 'LTP. Linearly towards inf.', curReds(2, :)};
methods = cell2struct(methods, {'name', 'title', 'color'}, 2);

info = Util.runInfo();
Util.showRunInfo(info);

%% Plot
rng(4);

areas = [1, 2, 3, numel(methods)];
areas = 10 .^ normrnd(ccMeanLog, ccStdLog, areas);
areas = reshape(num2cell(areas, 2), [], numel(methods));

curFig = figure();
curAx = axes(curFig);
hold(curAx, 'on');

for curAreaIdx = 1:size(areas, 1)
    for curMethodIdx = 1:numel(methods)
        curMethod = methods(curMethodIdx).name;
        curColor = methods(curMethodIdx).color;
    
        curAreas = ...
            connectEM.Consistency.Simulation.evolution( ...
                areas{curAreaIdx, curMethodIdx}, 'method', curMethod);
        areas{curAreaIdx, curMethodIdx} = curAreas;

        curX = std(curAreas, 0, 2) ./ mean(curAreas, 2);
        curY = log10(mean(curAreas, 2));
        
        scatter( ...
            curAx, curX(1), curY(1), 'o', ...
            'MarkerEdgeColor', 'none', ...
            'MarkerFaceColor', curColor);
        plot(curAx, curX, curY, 'Color', curColor);
    end
end

xlim(curAx, limX);
ylim(curAx, limY);

curLines = findobj(curAx, 'type', 'line');
set(curLines, 'LineWidth', 2);

xlabel(curAx, 'Coefficient of variation');
ylabel(curAx, 'log_{10}(Average ASI area [µm²])');
axis(curAx, 'square');

curLeg = flip(curLines(1:size(methods, 1)));
curLeg = legend(curLeg, {methods.title}, 'Location', 'SouthOutside');
curLeg.Box = 'off';

curFig.Position(3:4) = [315, 390];
connectEM.Figure.config(curFig, info);

%% Plot synapse pair trajectories
clear cur*;

curLimX = [1, 125];
curLimY = [0, 1];

curFig = figure();

for curAreaIdx = 1:size(areas, 1)
    for curMethodIdx = 1:numel(methods)
        curAreas = areas{curAreaIdx, curMethodIdx};
        
        curPlotIdx = ...
            curMethodIdx + ...
            numel(methods) * (curAreaIdx - 1);
        
        subplot(size(areas, 1), numel(methods), curPlotIdx);
        plot(curAreas);
    end
end


curAxes = curFig.Children;
set(curAxes, ...
    'XLim', curLimX, 'XTick', [], ...
    'YLim', curLimY, 'YTick', []);

curAx = curAxes(numel(methods));
xlabel(curAx, 'Time');
ylabel(curAx, 'Areas');

curLines = findobj(curFig, 'Type', 'Line');
set(curLines, 'LineWidth', 2);

connectEM.Figure.config(curFig);
curFig.Position(3:4) = [800, 160];

%% Plot average size and relative difference trajectory
clear cur*;

curLimX = [1, 125];
curLimLeftY = [0, 1];
curLimRightY = [0, 2];

curFig = figure();

for curAreaIdx = 1:size(areas, 1)
    for curMethodIdx = 1:numel(methods)
        curAreas = areas{curAreaIdx, curMethodIdx};
        
        curAvgArea = mean(curAreas, 2);
        curRelDiff = abs(diff(curAreas, 1, 2)) ./ curAvgArea;
        
        curPlotIdx = ...
            curMethodIdx + ...
            numel(methods) * (curAreaIdx - 1);
        
        curAx = subplot(size(areas, 1), numel(methods), curPlotIdx);
        yyaxis(curAx, 'left'); plot(curAvgArea);
        yyaxis(curAx, 'right'); plot(curRelDiff, 'k');
    end
end


curAxes = curFig.Children;
set(curAxes, 'XLim', curLimX, 'XTick', []);
arrayfun(@(ax) yyaxis(ax, 'right'), curAxes);
set(curAxes, 'YLim', curLimRightY, 'YTick', []);
arrayfun(@(ax) yyaxis(ax, 'left'), curAxes);
set(curAxes, 'YLim', curLimLeftY, 'YTick', []);

curAx = curAxes(numel(methods));
xlabel(curAx, 'Time');
ylabel(curAx, 'A / RD');

curLines = findobj(curFig, 'Type', 'Line');
set(curLines, 'LineWidth', 2);

connectEM.Figure.config(curFig);
curFig.Position(3:4) = [800, 160];
