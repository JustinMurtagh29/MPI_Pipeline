%% PROBLEM INSPECTOR
id = '20130916T234116';
load(['D:\sync\problemInspector\' id '.mat']);
vShort = v;
for vp=1:length(vShort)
    for i=1:length(vShort(vp).missions)
        firstFrames = cellfun(@length, {vShort(vp).missions(i).possibleEnds(:).firstFrame});
        lastFrames = cellfun(@length, {vShort(vp).missions(i).possibleEnds(:).lastFrame});
        keep = firstFrames | lastFrames;
        vShort(vp).missions(i).possibleEnds = vShort(vp).missions(i).possibleEnds(keep);
    end
end
colors = [1 0 0; 0 1 0; distinguishable_colors(max(cellfun(@length, {v(1).missions(:).possibleEnds}))-1, [1 0 0; 0 1 0])];

%% GLOBAL STATISTICS
visualizeRenderingLength(vShort, id);

%% Mission Visualization, Statistics, Videos similar to levelcreator
for i=1:50
    visualizeTasksNew(seg, vShort, parameter, kdbParameter, id, i, colors);
    arbitraryReslice(raw, seg, mito, vShort, parameter, kdbParameter, id, i, colors);
end

%% Stuff for Thomas/Levelcreator check

segNew = zeros(size(seg));
segNew(seg == missions(1).start.id) = 1;
segNew(seg == missions(1).end.id) = 2;
for i=1:length(missions(1).possibleEnds)
    if missions(1).possibleEnds(i).id ~= missions(1).end.id
        segNew(seg == missions(1).possibleEnds(i).id) = 1 + i;
    end
end

%%
seg = imdilate(seg, ones(3,3,3));
errorCenter = missions(1).errorCenter - parameter.bboxBig(:,1)';
direction = missions(1).start.direction;

%%
[a,b] = interpolateThis(raw, segNew, errorCenter, direction, [300 300], [50 50], [11.28 11.28 28], [11.28 11.28 11.28]);

%%
makeTaskMovieSetColors(uint16(b), a, ['C:\Users\mberning\Desktop\problemInspector\' id 'o' num2str(1, '%.3i') '.avi'], colors);
