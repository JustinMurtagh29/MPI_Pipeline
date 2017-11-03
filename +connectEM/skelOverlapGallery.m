function ovSkels = skelOverlapGallery( outputFolder, skels, agglos, ...
    skelToAgglos, ovT, segmentCom, noPlot, skelL, ...
    ovSkelLength, numOvR, erl, wavgL, recL )
%SKELOVERLAPGALLERY Produce a gallery of skeletons and their overlapping
%agglos.
% INPUT outputFolder: string
%           Save folder for the output files. If empty the figures are not
%           saved.
%       skels: [Nx1] cell of skeleton objects
%           The skeleton objects.
%       agglos: [Nx1] struct or [Nx1] cell
%           The agglos in the agglo or superagglo format. If the agglo
%           format is supplied then the segmentCom input must also be
%           supplied.
%       skelToAgglos: [Nx1] cell
%           Skeleton agglo overlap.
%           (see first output of L4.Agglo.aggloSkelOverlap)
%       ovT: (Optional) int
%           Minimal number of skeleton nodes in an agglo to consider the
%           agglo overlapping.
%           (Default: 10)
%       segmentCom: (Optional) [Nx3] int
%           The segment center-of-mass coordinates. This is only used if
%           the sagglos input is in the agglo format.
%       noPlot: (Optional) logical
%           Only return the outputs and do not plot the skeletons.
%           (Default: false)
%       skelL: (Optional) [Nx1] int
%           Length of the ground truth skeletons for the output in um.
%           (Default: will be calculated)
%       ovSkelLength: (Optional) [Nx1] cell of [Mx1] double
%           The length of the agglos that is actually overlapping with the
%           ground truth skeleton.
%       numOvR: (Optional) [Nx1] int
%           The number of segments/agglos overlapping with the
%           corresponding skeleton that are not dendrite agglos.
%       erl: (Optional) [Nx1] double
%           The erl in um for the corresponding skeleton.
%       wavgL: (Optional) [Nx1] double
%           The weighted average agglo length in um for the corresponding
%           skeleton.
%       recL: (Optional) [Nx1] double
%           The total recovered ground truth path length.
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

if ~isempty(outputFolder) && ~exist(outputFolder, 'dir')
    mkdir(outputFolder)
end
if ~isempty(outputFolder)
    Util.log('Gallery is stored at %s.', outputFolder);
end

if ~exist('ovT', 'var') || isempty(ovT) || ~isnumeric(ovT)
    ovT = 10;
end
ovT = round(ovT);
Util.log('ovT is set to %d.', ovT);

if ~exist('noPlot', 'var') || isempty(noPlot)
    noPlot = false;
end

if ~exist('skelL', 'var') || isempty(skelL)
    skelL = cellfun(@(x)x.pathLength()./1000, skels);
end
if ~exist('ovSkelLength', 'var') || isempty(ovSkelLength)
    ovSkelLength = cellfun(@(x)nan(size(x, 1), 1), skelToAgglos, 'uni', 0);
end
if ~exist('numOvR', 'var') || isempty(numOvR)
    numOvR = NaN(length(skels), 1);
end
if ~exist('erl', 'var') || isempty(erl)
    erl = NaN(length(skels), 1);
end
if ~exist('wavgL', 'var') || isempty(wavgL)
    wavgL = NaN(length(skels), 1);
end
if ~exist('recL', 'var') || isempty(recL)
    recL = NaN(length(skels), 1);
end

% add overlapping agglos to gt skeletons
ovSkels = cell(length(skels), 1);
numGTTrees = cellfun(@(x)x.numTrees(), skels);
N_ov = zeros(length(skels), 1);
for i = 1:length(skels)
    skel = skels{i};
    skel.verbose = false;
    if numGTTrees(i) == 1
        skel.names{1} = 'GroundTruthSkeleton';
    else
        skel.names = arrayfun(@(x) ...
            sprintf('GroundTruthSkeleton_Part%d', x), numGTTrees, ...
            'uni', 0);
    end
    skel = skel.setDescription(sprintf(['Whole cell ground truth ' ...
        'tracing (%s) with overlapping agglos (ovT %d).'], ...
        skel.filename, ovT));
    ov = skelToAgglos{i}(:,2);
    idx = skelToAgglos{i}(ov >= ovT, 1);
    N_ov(i) = length(idx);
    if isstruct(agglos)
        skel = L4.Agglo.superAgglo2Skel(agglos(idx), skel);
    else
        skel = L4.Agglo.agglo2Nml(agglos(idx), segmentCom, skel);
    end
    skel.names(end - length(idx) + 1:end) = arrayfun(@(x, y) ...
            sprintf('Agglo_%d_ov%d', x, y), idx, ov(ov >= ovT), ...
            'uni', 0);
    ovSkels{i} = skel;
end

if noPlot
    return;
end

% plotting
bbox = [129 5574; 129 8509; 129 3414];
bbox = bsxfun(@times, bbox, [11.24; 11.24; 28]./1000);
for i = 1:length(ovSkels)
    %set f to din A4 size
    f = figure;
    f.Units = 'centimeters';
    f.Position = [1 1 21 29.7];
    Visualization.Figure.adaptPaperPosition();
    for j = 1:4
        a = subplot(2, 2, j);
        ovSkels{i}.plot(1:numGTTrees(i), [0, 0, 0], true, 2);
        hold on
        ovSkels{i}.plot(numGTTrees(i)+1:ovSkels{i}.numTrees(), hsv(6), true);
        xlabel('x (um)');
        ylabel('y (um)');
        zlabel('z (um)');
        grid on
        switch j
            case 1
                %xy
                view([0, 90]);
                filename = strrep(ovSkels{i}.filename, '_', '\_');
                txt = sprintf(['%s\n' ...
                    'GT Skel length: %.2f um; recovered Length: %.2f um; recall: %.2f\n' ...
                    'Number of overlaps: %d (dendrites) %d (rest)\n', ...
                    'Overlap length: %s um\n', ...
                    'ERL: %.1f um; WeighAvgL: %.1f um'], ...
                    filename, skelL(i), recL(i), recL(i) / skelL(i), ...
                    N_ov(i), numOvR(i), ...
                    mat2str(round(100.*ovSkelLength{i})/100), ...
                    erl(i), wavgL(i));

                text(a, 0, 1, txt, ...
                    'FontSize', 8, 'FontName', 'Arial', ...
                    'HorizontalAlignment', 'left', ...
                    'Units', 'normalized', ...
                    'VerticalAlignment', 'bottom');
            case 2
                %xz
                view([0, 0]);
            case 3
                %yz
                view([90, 0]);
            case 4
                % default 3d view
        end
        axis equal
        a.XLim = bbox(1,:);
        a.YLim = bbox(2,:);
        a.ZLim = bbox(3,:);
    end
    if ~isempty(outputFolder)
        if isempty(ovSkels{i}.filename)
            ovSkels{i}.filename = sprintf('Skel%03d', i);
        end
        outfile = fullfile(outputFolder, ...
            sprintf('Skel%02d_%s_AggloOverlap.pdf', i, ovSkels{i}.filename));
        print(outfile, '-dpdf')
        Util.log('Saving file %s.', outfile);
        close(f);
    end
end

end
