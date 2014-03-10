%% Load data
clear all; clc;
cubes = [13 15; 14 16; 17 19];
coords(:,1) = 1 + 128 * cubes(:,1);
coords(:,2) = 128 + 128 * cubes(:,2);
raw = readKnossosRoi('E:\e_k0563\k0563_mag1\', '100527_k0563_mag1', coords);
load('P:\segForPaper.mat');
skel = readNml('I:\CortexConnectomics\Manuel\backup\20130726Laptop\c\data\ssdf.221.nmllocal');
seg = permute(seg, [2 1 3]);

%% Calculate isosurfaces of all supervoxels
issfs = cell(max(seg(:)),1);
seg = padarray(seg, [3 3 3]);
for i=1:max(seg(:))
    issfs{i} = isosurface(smooth3(seg == i, 'gaussian', [5 5 5], 1), .9);
    if ~isempty(issfs{i}.vertices)
        issfs{i}.vertices = bsxfun(@minus, issfs{i}.vertices, [3 3 3]);  
    end
    display(num2str(i, '%.4i'));
end
seg = seg(4:end-3,4:end-3,4:end-3);
save('C:\Users\mberning\Desktop\201310017issfsForPaperFigure6.mat');

%% Find relation between skeletons and segments (necessary if coloring segments according to skeleton, probably not a good idea)
% equivMatrix = zeros(size(skel,2), single(max(seg(:))));
% for l=1:length(skel)
%     nodes{l} = skel{l}.nodes(:,1:3);
%     for m=1:size(nodes{l},1)
%         if seg(nodes{l}(m,1), nodes{l}(m,2), nodes{l}(m,3))
%             equivMatrix(l,seg(nodes{l}(m,1), nodes{l}(m,2), nodes{l}(m,3))) = equivMatrix(l,seg(nodes{l}(m,1), nodes{l}(m,2), nodes{l}(m,3))) + 1;
%         end
%     end
% end
% equivMatrixBinary = equivMatrix >= 1;
% vec = sum(equivMatrixBinary,2);
% idx = find(vec > 1);
% obj = cell(length(idx),1);
% for m=1:length(idx)
%     obj{m} = find(equivMatrixBinary(idx(m),:));
% end

%% Some settings
colors = distinguishable_colors(max(seg(:)), [0 0 0]);

%% Figure segmentation
close all;
figure('Units', 'centimeters', 'Position', [0 0 29.7 21], 'Renderer', 'OpenGL', ...
    'PaperType', 'A4', 'PaperUnits', 'centimeters', 'PaperPosition', [0 0 29.7 21], 'PaperOrientation', 'portrait');
hold on;
for i=1:length(issfs)
    hpatch = patch(issfs{i});
    set(hpatch,'FaceColor', colors(i,:), 'EdgeColor', 'none');
    set(hpatch,'AmbientStrength', .5, 'SpecularStrength', .5, 'DiffuseStrength', .5);
end

slices = {[384] [384] [1]};
plotOriginalData(double(raw), slices);

caxis([50 200]);
colormap('gray');

lighting gouraud;
camproj perspective;
daspect([25 25 12]);
view(3);
hlight = camlight('headlight');
axis off;

set(gcf,'PaperSize',fliplr(get(gcf,'PaperSize')));
print(gcf, 'C:\Users\mberning\Desktop\figure6_1.png', '-dpng', '-r900');

%% Figure skeletons
figure('Units', 'centimeters', 'Position', [0 0 29.7 21], 'Renderer', 'OpenGL', ...
    'PaperType', 'A4', 'PaperUnits', 'centimeters', 'PaperPosition', [0 0 29.7 21], 'PaperOrientation', 'portrait');
plotSkelTraceNips(raw, skel, colors);

set(gcf,'PaperSize',fliplr(get(gcf,'PaperSize')));
print(gcf, 'C:\Users\mberning\Desktop\figure6_2.png', '-dpng', '-r900');
