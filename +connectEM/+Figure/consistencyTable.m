% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

info = Util.runInfo();
Util.showRunInfo(info);

%% Configuration
rowNames = {'Exc', 'CC', 'TC'};
colNames = {'All', 'WC', 'AD', 'OD'};
showNumbers = true;
colorN = 256;

% First dimension corresponds to small (1) and large (2) low-CV regions
data = nan(2, numel(rowNames), numel(colNames));

% Numbers from connectEM.Connectome plotSynapseSizeConsistency
% Commit 078e1a081c259ffb956144dfc3ddc595710e02e2
data(:, 1, 1) = [15.83, 16.23];
data(:, 1, 2) = [13.44, 12.11];
data(:, 1, 3) = [ 5.69,  8.77];
data(:, 1, 4) = [12.29, 13.17];

data(:, 2, 1) = [15.64, 17.42];
data(:, 2, 2) = [13.14, 11.73];
data(:, 2, 3) = [  nan,   nan];
data(:, 2, 4) = [13.47, 11.35];

data(:, 3, 1) = [ 2.05,  6.61];
data(:, 3, 2) = [  nan,   nan];
data(:, 3, 3) = [  nan,   nan];
data(:, 3, 4) = [ 0.66,  3.72];

colors = reshape(linspace(0, 1, colorN), [], 1);
colors = { ...
    (ones(1, 3) - colors) + colors .* [0.4660, 0.6740, 0.1880]; ...
    (ones(1, 3) - colors) + colors .* [0.8500, 0.3250, 0.0980]};
assert(isequal(numel(colors), size(data, 1)));

%% Plot
curCaxisMax = max(data(:), [], 'omitnan');
curCaxisMax = ceil(curCaxisMax / 10) * 10;

close all;
fig = figure();
ax = axes(fig);
hold(ax, 'on');

for curRow = 1:numel(rowNames)
    curRowName = rowNames{curRow};
    
    for curCol = 1:numel(colNames)
        curColName = colNames{curCol};
        
        curCoord = nan(3, 2);
        curCoord(:, 1) = (curCol - 1) + [0, 0, 1] + 0.5;
        curCoord(:, 2) = (curRow - 1) + [0, 0, 1] + 0.5;
        
        for curId = 1:2
            curData = data(curId, curRow, curCol);
            if isnan(curData); continue; end
            
            switch curId
                case 1, curCoord(2, :) = [curCoord(1, 1), curCoord(3, 2)];
                case 2, curCoord(2, :) = [curCoord(3, 1), curCoord(1, 2)];
                otherwise, error('Invalid id "%d"', curId);
            end
            
            curShape = polyshape(curCoord(:, 1), curCoord(:, 2));
            
            curColor = linspace(0, curCaxisMax, colorN);
            curColor = colors{curId}(discretize(curData, curColor), :);
            
            plot(curShape, 'FaceColor', curColor, 'FaceAlpha', 1);
            
            text( ...
                ax, ...
                mean(curCoord(:, 1)), ...
                mean(curCoord(:, 2)), ...
                sprintf('%.0f', curData), ...
                'HorizontalAlignment', 'center');
        end
    end
end

axis('equal');
ax.YDir = 'reverse';
ax.XAxisLocation = 'top';

xlim(ax, [0, numel(colNames)] + 0.5);
xticks(ax, 1:numel(colNames));
xticklabels(ax, colNames);

ylim(ax, [0, numel(rowNames)] + 0.5);
yticks(ax, 1:numel(rowNames));
yticklabels(ax, rowNames);

fig.Position(3:4) = [310, 250];
connectEM.Figure.config(fig, info);
