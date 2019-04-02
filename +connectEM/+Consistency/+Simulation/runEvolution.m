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

% Run simulations
rng(4);

areas = [1, 2, 3, numel(methods)];
areas = 10 .^ normrnd(ccMeanLog, ccStdLog, areas);
areas = reshape(num2cell(areas, 2), [], numel(methods));

for curAreaIdx = 1:size(areas, 1)
    for curMethodIdx = 1:numel(methods)
        curMethod = methods(curMethodIdx).name;
    
        areas{curAreaIdx, curMethodIdx} = ...
            connectEM.Consistency.Simulation.evolution( ...
                areas{curAreaIdx, curMethodIdx}, 'method', curMethod);
    end
end

%% Plot
rng(4);

curFig = figure();
curAx = axes(curFig);
hold(curAx, 'on');

for curAreaIdx = 1:size(areas, 1)
    for curMethodIdx = 1:numel(methods)
        curMethod = methods(curMethodIdx).name;
        curColor = methods(curMethodIdx).color;
        curAreas = areas{curAreaIdx, curMethodIdx};

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

curLimX = [0, 240];
curTicksX = linspace(curLimX(1), curLimX(2), 4);
curLimY = [0, 1.5];
curTicksY = [0, 1];

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
    'XLim', curLimX, 'XTick', curTicksX, ...
    'YLim', curLimY, 'YTick', curTicksY);

curTickLabels = xticklabels(curAxes(1));
curTickLabels(2:(end - 1)) = {''};
xticklabels(curAxes, curTickLabels);

curAx = curAxes(numel(methods));
xlabel(curAx, 'Time step');
ylabel(curAx, 'Area [µm²]');

curLines = findobj(curFig, 'Type', 'Line');
set(curLines, 'LineWidth', 2);

connectEM.Figure.config(curFig);
curFig.Position(3:4) = [315, 215];

%% Movie
clear cur*;

curSteps = size(areas{1}, 1);
curFrames = struct('cdata', {}, 'colormap', {});

curFig = figure();
curAx = axes(curFig);
hold(curAx, 'on');

xlim(curAx, [0, 2]);
ylim(curAx, limY);

axis(curAx, 'square');
xticklabels(curAx, {});
yticklabels(curAx, {});

curFig.Position(3:4) = 200;
connectEM.Figure.config(curFig);

for curT = 1:size(areas{1}, 1)
    delete(curAx.Children);
    
    for curAreaIdx = 1:size(areas, 1)
        for curMethodIdx = 1:numel(methods)
            curColor = methods(curMethodIdx).color;
            curAreas = areas{curAreaIdx, curMethodIdx};

            curX = abs(diff(curAreas, 1, 2)) ./ mean(curAreas, 2);
            curY = log10(mean(curAreas, 2));

            scatter( ...
                curAx, curX(1), curY(1), 'o', ...
                'MarkerEdgeColor', 'none', ...
                'MarkerFaceColor', curColor);
            plot(curAx, ...
                curX(1:curT), curY(1:curT), ...
                'Color', curColor, ...
                'LineWidth', 2);
        end
    end
    
    curFrames(curT) = getframe(curFig);
end
