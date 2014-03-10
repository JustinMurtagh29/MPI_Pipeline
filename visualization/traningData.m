% Visualize training data
%clear all; clc;
files = dir('D:\sync\trainingData\2013*.mat');
load D:\sync\trainingData\exclude.mat;
load D:\sync\trainingData\parameter.mat;
for i=1:length(files);
    if any(stacks(i).taskID == excludeTask)
        exclude(i) = 1;
    else
        exclude(i) = 0;
    end
end




%% Find correspondence of index shifts due to exclusion of stacks
indShifts = cumsum(~exclude);
find(indShifts == 67)

%% 
for i=1:219
    load(['D:\sync\trainingData\' files(i).name]);
    colors = distinguishable_colors(max(stack(:)), [0 0 0]);
    stack = permute(stack, [2 1 3]);
    figure('Units', 'centimeters', 'Position', [0 0 29.7 21], 'Renderer', 'OpenGL', ...
    'PaperType', 'A4', 'PaperUnits', 'centimeters', 'PaperPosition', [0 0 29.7 21], 'PaperOrientation', 'portrait');
    hold on;
    uniqueValues = unique(stack);
    uniqueValues(uniqueValues == 0) = [];
    for j=1:length(uniqueValues)
        issf = isosurface(smooth3(padarray(stack == uniqueValues(j), [3 3 3], 0), 'gaussian', [5 5 5], 1.5), .2);
        if length(issf.vertices) > 1
            issf.vertices = bsxfun(@minus, issf.vertices, [3 3 3]);
            hpatch = patch(issf);
            set(hpatch,'FaceColor', colors(uniqueValues(j),:), 'EdgeColor', 'none');
            set(hpatch,'AmbientStrength', .2, 'SpecularStrength',.6, 'DiffuseStrength',.6);
        end
    end
    [x, y, z] = meshgrid(1, [1 100], [1 100]);
    x = squeeze(x); y = squeeze(y); z = squeeze(z); 
    surf(x, y, z, squeeze(raw(1,1:100,:)), 'Facecolor', 'texture', 'EdgeColor', 'none');
    [x, y, z] = meshgrid([1 100], 1, [1 100]);
    x = squeeze(x); y = squeeze(y); z = squeeze(z); 
    surf(x, y, z, squeeze(raw(1:100,1,:)), 'Facecolor', 'texture', 'EdgeColor', 'none');
    [x, y, z] = meshgrid([1 100], [1 100], 1);
    x = squeeze(x); y = squeeze(y); z = squeeze(z); 
    surf(x, y, z, squeeze(raw(1:100,1:100,1)'), 'Facecolor', 'texture', 'EdgeColor', 'none');
    caxis([50 200]);
    colormap('gray');
    lighting gouraud;
    camproj perspective;
    view(113,14);
    hlight = camlight('headlight');
    daspect([28 28 11.24]);
    axis off;
    stack = permute(stack, [2 1 3]);
    if exclude(i)
        saveas(gcf, ['D:\sync\trainingData\issfs' num2str(i,'%.3i') 'EXC.fig']);
        makeTaskMovieSetColors(stack, raw, ['D:\sync\trainingData\seg' num2str(i,'%.3i') 'EXC'], colors);
        makeTargetMovie(target, raw, ['D:\sync\trainingData\target' num2str(i,'%.3i') 'EXC']);   
    else
        saveas(gcf, ['D:\sync\trainingData\issfs' num2str(i,'%.3i') '.fig']);
        makeTaskMovieSetColors(stack, raw, ['D:\sync\trainingData\seg' num2str(i,'%.3i')], colors);
        makeTargetMovie(target, raw, ['D:\sync\trainingData\target' num2str(i,'%.3i')]);
    end
    close all;
end

%% For paper Figure 2
% Debug point in line 6- and go there
f = 91;
plot3([1 100], [100 100], [f f], 'c', 'LineWidth', 2);
plot3([1 100], [1 1], [f f], 'c', 'LineWidth', 2);
plot3([1 1], [1 100], [f f], 'c', 'LineWidth', 2);
plot3([100 100], [1 100], [f f], 'c', 'LineWidth', 2);

figure;
img = stack(:,:,f);
imshow(raw(:,:,f), [40 210]);
hold on;
temp = label2rgb(img, colors, [0 0 0]);
himage = imshow(temp);
set(himage, 'AlphaData', 0.4);

%% 
addpath(genpath('C:/code/KLEE/'));
i = 187;
load(['D:\sync\trainingData\' files(i).name]);
KLEE_v4('stack', raw, 'stack_2', target);

%%
for i=1:279
    load(['D:\sync\trainingData\' files(i).name]);
    voxelIntra(i) = sum(target(:) == 1);
    voxelExtra(i) = sum(target(:) == -1);
    voxelUnlabeled(i) = sum(target(:) == 0);
    meanValue(i) = mean(raw(:));
    raw = single(raw) - 122;
    stdValue(i) = std(raw(:));
%     pause(.2);
end

%%
figure;
plot([voxelIntra; voxelExtra; voxelUnlabeled]');
hold on;
for i=1:length(exclude)
    if exclude(i)
        plot([i i], [0 1e6], 'c');
    end
end
legend('Intracellular voxel', 'Extracellular voxel', 'Unlabeled voxel', 'Excluded Tasks');

%%
figure;
plot([meanValue; stdValue]');
hold on;
plot([0 279], [122 122], 'r');
plot([0 279], [22 22], 'y');
for i=1:length(exclude)
    if exclude(i)
        plot([i i], [0 140], 'c');
    end
end
legend('Mean Raw Values', 'Standard deviation raw values', 'Mean Normalization', 'Standard Deviation Normalization', 'Excluded Tasks');
