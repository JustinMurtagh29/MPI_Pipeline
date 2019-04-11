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
    'ltdLinearZero', 'LTD. Linearly towards zero.', nan(1, 3); ...
    'ltdLinearNonZero', 'LTD. Linearly towards non-zero.', curReds(1, :);
    'ltdExponentialZero', 'LTD. Exponential decay to zero.', nan(1, 3); ...
    'ltdExponentialNonZero', 'LTD. Exponential decay to non-zero.', curReds(2, :); ...
    'ltpExponentialFinite', 'LTP. Exponential grows to maximum.', curGreens(2, :); ...
    'ltpLinearInf', 'LTP. Linearly towards inf.', curGreens(4, :)};
methods = cell2struct(methods, {'name', 'title', 'color'}, 2);

info = Util.runInfo();
Util.showRunInfo(info);

% Run simulations
rng(4);

areas = [1, 2, 3, numel(methods)];
areas = 10 .^ normrnd(ccMeanLog, ccStdLog, areas);
areas = reshape(num2cell(areas, 2), [], numel(methods));

for curAreaIdx = 1:size(areas, 1)
    for curMethodId = 1:numel(methods)
        curMethod = methods(curMethodId).name;
    
        areas{curAreaIdx, curMethodId} = ...
            connectEM.Consistency.Simulation.evolution( ...
                areas{curAreaIdx, curMethodId}, 'method', curMethod);
    end
end

%% Plot
rng(4);

curFig = figure();
curAx = axes(curFig);
hold(curAx, 'on');

for curAreaIdx = 1:size(areas, 1)
    for curMethodId = 1:numel(methods)
        curMethod = methods(curMethodId).name;
        curColor = methods(curMethodId).color;
        curAreas = areas{curAreaIdx, curMethodId};

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
    for curMethodId = 1:numel(methods)
        curAreas = areas{curAreaIdx, curMethodId};
        
        curPlotIdx = ...
            curMethodId + ...
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

%% Movie of synapse pair trajectories
clear cur*;

curOutFile = '/home/amotta/Desktop/size-vs-time_movie-%d';

curPairConfigs = @(r, c) sub2ind(size(areas), r, c);
curPairConfigs = { ...
    [curPairConfigs(2, 5), curPairConfigs(3, 6)], ...
    [curPairConfigs(1, 4), curPairConfigs(1, 2)]};

for curPairConfigIdx = 1:numel(curPairConfigs)
    curPairIds = curPairConfigs{curPairConfigIdx};
    curFrames = struct('cdata', {}, 'colormap', {});
    
    curFig = figure();
    curFig.Position(3:4) = [430, 240];
    
    for curT = 1:size(areas{1}, 1)
        delete(curFig.Children);
        
        for curPairIdx = 1:numel(curPairIds)
            curPairId = curPairIds(curPairIdx);
            curAreas = areas{curPairId};
            
            curAx = subplot(1, numel(curPairIds), curPairIdx);
            hold(curAx, 'on');
            
            plot(curAx, ...
                1:curT, transpose(curAreas(1:curT, 1)), ...
                'Color', 'black', 'LineWidth', 2);
            plot(curAx, ...
                1:curT, transpose(curAreas(1:curT, 2)), ...
                'Color', 'black', 'LineWidth', 2);
            
            xlim(curAx, [0, size(areas{1}, 1)]);
            ylim(curAx, [0, 1]);
            
            xlabel(curAx, []);
            xticks(curAx, []);
            
            ylabel(curAx, []);
            yticks(curAx, [0, 1]);
            yticklabels(curAx, {});
        end
        
        connectEM.Figure.config(curFig);
        curFrames(end + 1) = getframe(curFig); %#ok
    end
    
    curWriter = sprintf(curOutFile, curPairConfigIdx);
    curWriter = VideoWriter(curWriter, 'Uncompressed AVI'); %#ok
    curWriter.open();
    curWriter.writeVideo(curFrames);
    curWriter.close();
end

%% Movie
clear cur*;

curOutFile = '/home/amotta/Desktop/cv-vs-log10-avg-size_movie-%d';
curMethodConfigs = {[5, 6], [5, 6, 2, 4]};

for curMethodConfigIdx = 1:numel(curMethodConfigs)
    curMethodIds = curMethodConfigs{curMethodConfigIdx};
    curFrames = struct('cdata', {}, 'colormap', {});

    curFig = figure();
    curAx = axes(curFig); %#ok
    hold(curAx, 'on');

    xlim(curAx, [0, 1.5]);
    ylim(curAx, limY);

    axis(curAx, 'square');
    xticklabels(curAx, {});
    yticklabels(curAx, {});

    curFig.Position(3:4) = 600;
    connectEM.Figure.config(curFig);

    for curT = 1:2:size(areas{1}, 1)
        delete(curAx.Children);

        for curAreaIdx = 1:size(areas, 1)
            for curMethodId = curMethodIds
                curColor = methods(curMethodId).color;
                curAreas = areas{curAreaIdx, curMethodId};
                
                curX = std(curAreas, 0, 2) ./ mean(curAreas, 2);
                curY = log10(mean(curAreas, 2));

                scatter( ...
                    curAx, curX(1), curY(1), 96, 'o', ...
                    'MarkerEdgeColor', 'none', ...
                    'MarkerFaceColor', curColor);
                plot(curAx, ...
                    curX(1:curT), curY(1:curT), ...
                    'Color', curColor, ...
                    'LineWidth', 3);
            end
        end
    
        curAx.XAxis.TickLength(1) = 0.02;
        curAx.XAxis.LineWidth = 2;
        curAx.YAxis.TickLength(1) = 0.02;
        curAx.YAxis.LineWidth = 2;

        curFrames(end + 1) = getframe(curFig);
    end
    
    curWriter = sprintf(curOutFile, curMethodConfigIdx);
    curWriter = VideoWriter(curWriter, 'Uncompressed AVI'); %#ok
    curWriter.open();
    curWriter.writeVideo(curFrames);
    curWriter.close();
end
