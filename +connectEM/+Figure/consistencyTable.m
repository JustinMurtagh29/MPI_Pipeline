% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

info = Util.runInfo();
Util.showRunInfo(info);

%% Configuration
rowNames = {'Exc', 'CC', 'TC'};
colNames = {'All', 'PD', 'AD', 'OD'};
colorN = 256;

% First dimension corresponds to minimum and maximum.
% Second dimension corresponds to small (1) and large (2) low-CV regions.
data = nan(2, 2, numel(rowNames), numel(colNames));

% Numbers from connectEM.Connectome plotSynapseSizeConsistency
% Commit 2574326f1e904c8c68a080a60dfc43b567e9b9e1

% Excitatory synapses onto all targets
data(:, 1, 1, 1) = [15.2, 18.6];
data(:, 2, 1, 1) = [15.8, 19.7];
% Excitatory synapses onto proximal dendrites
data(:, 1, 1, 2) = [12.8, 17.9];
data(:, 2, 1, 2) = [11.4, 17.9];
% Excitatory synapses onto apical dendrites
data(:, 1, 1, 3) = [ 5.1,  7.4];
data(:, 2, 1, 3) = [ 8.2,  9.7];
% Excitatory synapses onto other dendrites
data(:, 1, 1, 4) = [11.5, 17.9];
data(:, 2, 1, 4) = [12.6, 19.9];

% Corticocortical synapses onto all targets
data(:, 1, 2, 1) = [15.2, 20.7];
data(:, 2, 2, 1) = [16.9, 22.0];
% Corticocortical synapses onto proximal dendrites
data(:, 1, 2, 2) = [12.7, 17.6];
data(:, 2, 2, 2) = [11.5, 17.5];
% Corticocortical synapses onto apical dendrites
data(:, :, 2, 3) = nan;
% Corticocortical synapses onto other dendrites
data(:, 1, 2, 4) = [12.8, 18.5];
data(:, 2, 2, 4) = [10.5, 19.7];

% Thalamocortical synapses onto all targets
data(:, 1, 3, 1) = [ 1.7,  6.7];
data(:, 2, 3, 1) = [ 6.2, 15.4];
% Thalamocortical synapses onto proximal dendrites
data(:, :, 3, 2) = nan;
% Thalamocortical synapses onto apical dendrites
data(:, :, 3, 3) = nan;
% Thalamocortical synapses onto other dendrites
data(:, 1, 3, 4) = [ 0.3,  5.8];
data(:, 2, 3, 4) = [ 3.4, 16.7];

colors = reshape(linspace(0, 1, colorN), [], 1);
colors = { ...
    (ones(1, 3) - colors) + colors .* [0.4660, 0.6740, 0.1880]; ...
    (ones(1, 3) - colors) + colors .* [0.8500, 0.3250, 0.0980]};
assert(isequal(numel(colors), size(data, 1)));

%% Plot
curCaxisMax = mean(reshape(data, 2, []));
curCaxisMax = max(curCaxisMax, [], 'omitnan');
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
            curData = data(:, curId, curRow, curCol);
            if all(isnan(curData)); continue; end
            curMean = mean(curData);
            
            switch curId
                case 1, curCoord(2, :) = [curCoord(1, 1), curCoord(3, 2)];
                case 2, curCoord(2, :) = [curCoord(3, 1), curCoord(1, 2)];
                otherwise, error('Invalid id "%d"', curId);
            end
            
            curShape = polyshape(curCoord(:, 1), curCoord(:, 2));
            
            curColor = linspace(0, curCaxisMax, colorN);
            curColor = discretize(mean(curData), curColor);
            curColor = colors{curId}(curColor, :);
            
            plot(curShape, 'FaceColor', curColor, 'FaceAlpha', 1);
            
            text( ...
                ax, ...
                mean(curCoord(:, 1)), ...
                mean(curCoord(:, 2)), ...
                sprintf('%.0f-%.0f', curData), ...
                'HorizontalAlignment', 'center');
        end
    end
end

colormap(ax, colors{2});
caxis(ax, [0, curCaxisMax]);
cbar = colorbar('peer', ax);
cbar.Label.String = 'Percent consistent';

axis('equal');
ax.YDir = 'reverse';
ax.XAxisLocation = 'top';

xlim(ax, [0, numel(colNames)] + 0.5);
xticks(ax, 1:numel(colNames));
xticklabels(ax, colNames);

ylim(ax, [0, numel(rowNames)] + 0.5);
yticks(ax, 1:numel(rowNames));
yticklabels(ax, rowNames);

fig.Position(3:4) = [420, 330];
connectEM.Figure.config(fig, info);
