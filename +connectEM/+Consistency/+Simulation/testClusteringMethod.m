% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
connFile = '/tmpscratch/amotta/l4/2019-02-15-consistency-analysis-using-gaussian-mixture-models/connection-types_v2.mat';
analysisFile = '/tmpscratch/amotta/l4/2019-02-15-consistency-analysis-using-gaussian-mixture-models/connection-types_weight-matrix-analysis_v1.mat';

info = Util.runInfo();
Util.showRunInfo(info);

%% Load data
load(connFile, 'asiT', 'preT', 'postT');
load(analysisFile, 'axonPerm', 'dendPerm');

[~, connName] = fileparts(connFile);

%% Plot precision / recall
clear cur*;

curPlots = struct;
curPlots(1).title = 'Axons';
curPlots(1).perm = axonPerm;
curPlots(1).classes = preT.classId(unique(asiT.preAggloId));

curPlots(2).title = 'Dendrite';
curPlots(2).perm = dendPerm;
curPlots(2).classes = postT.classId(unique(asiT.postAggloId));

curFig = figure();
curFig.Color = 'white';
curFig.Position(3:4) = [250, 225];

curAx = axes(curFig);
hold(curAx, 'on');

for curPlot = curPlots 
    assert(isequal(numel(curPlot.perm), numel(curPlot.classes)));
    
    curClasses = curPlot.classes(curPlot.perm);
    curCount = reshape(cumsum(curClasses == 1), 1, []);
    curPrec = curCount ./ (1:numel(curCount));
    curRec = curCount / curCount(end);
    
    plot(curAx, curRec, curPrec);
end

set(curAx.Children, 'LineWidth', 2);
curAx.TickDir = 'out';
axis(curAx, 'square');
grid(curAx, 'on');

xlim(curAx, [0, 1]);
xticks(curAx, linspace(0, 1, 6));
xlabel(curAx, 'Recall');

ylim(curAx, [0, 1]);
yticks(curAx, linspace(0, 1, 6));
ylabel(curAx, 'Precision');

curLeg = legend(curAx, {curPlots.title}, 'Location', 'SouthWest');
curLeg.Box = 'off';

title(curAx, ...
    {info.filename; info.git_repos{1}.hash; connName}, ...
    'FontWeight', 'normal', 'FontSize', 8, 'Interpreter', 'none');
