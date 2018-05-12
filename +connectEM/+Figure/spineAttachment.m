% The script to generate the isosurfaces in the first place is here:
%   +connectEM/+WholeCell/+Script/calculateIsosurfaces.m
%
% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
isoDir = '/tmpscratch/amotta/l4/2018-05-10-whole-cell-isosurfaces-spine-evolution';

inFile = fullfile(isoDir, '%2$s', 'mat', 'iso-%1$d.mat');
outFile = fullfile(isoDir, 'renderings', 'cell-%d_%s-spines_v1.png');

isoTags = {'pre', 'auto', 'full'};
isoCount = 96;

info = Util.runInfo();

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'), 'p');
param = param.p;

%% Process cells
for curCellId = 1:isoCount
    curIsos = cell(size(isoTags));
    for curIdx = 1:numel(isoTags)
        curIsoFile = sprintf(inFile, curCellId, isoTags{curIdx});

        curIso = load(curIsoFile);
        curIso = curIso.isoSurf;

        % To physical units. Make Z the cortical axis.
        curIso.vertices = curIso.vertices .* param.raw.voxelSize;
        curIso.vertices = curIso.vertices(:, [3, 2, 1]);
        
        curIsos{curIdx} = curIso;
    end

    renderCell(param, info, isoTags, outFile, curCellId, curIsos);
end

%% Utility
function renderCell(param, info, isoTags, outPath, cellId, isos)
    fig = figure();
    fig.Visible = 'off';
    fig.Color = 'black';
    fig.Position(3:4) = [2560, 1357];

    ax = axes(fig);
    ax.Color = 'black';
    axis(ax, 'equal');
    
    colors = ax.ColorOrder;
    colors = colors(1:numel(isoTags), :);
    colors = num2cell(colors, 2);

    patches = flip(cellfun( ...
        @(iso) patch(ax, iso), flip(isos)));
   [patches.FaceColor] = deal(colors{:});
   [patches.EdgeColor] = deal('none');
    material(patches, 'dull');

    % Fix Cortical axis
    ax.ZAxis.Direction = 'reverse';
    ax.YAxis.Direction = 'reverse';

    curAxes = [ax.XAxis, ax.YAxis, ax.ZAxis];
   [curAxes.Visible] = deal('off');
   [curAxes.Color] = deal('white');
   
    % Bounding box
    curLims = param.raw.voxelSize(:) .* param.bbox;
    curLims = num2cell(curLims([3, 2, 1], :), 2);
   [curAxes.Limits] = deal(curLims{:});

    ax.Box = 'on';
    ax.BoxStyle = 'back';
    
    % Show evolution
   [patches.Visible] = deal('off');
    patches(1).Visible = 'on';

    % View
    view(ax, 3);
    camlight(ax, -37.5, 30);
    
    for curId = 1:numel(isoTags)
       [patches(1:curId).Visible] = deal('on');
        title(ax, { ...
            info.filename; info.git_repos{1}.hash; ...
            sprintf('Cell %d (%s spines)', cellId, isoTags{curId})}, ...
            'Color', 'white', 'FontWeight', 'normal', 'FontSize', 10);
        curOutPath = sprintf(outPath, cellId, isoTags{curId});
        export_fig('-r172', curOutPath);
    end
    
    delete(fig);
end
