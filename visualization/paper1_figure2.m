load('I:\CortexConnectomics\Manuel\backup\20130726Laptop\c\data\toAmira_Neu.mat');
%%
for i = 1:4
    issfs{i} = isosurface(smooth3(trace == i, 'gaussian', [5 5 5], 1.5));
end
raw = permute(raw, [2 1 3]);

%%
close all;
colors = {'r', 'g', 'b', 'y'};
zPlane = 141;
figure('Units', 'centimeters', 'Position', [0 0 29.7 21], 'Renderer', 'OpenGL', ...
    'PaperType', 'A4', 'PaperUnits', 'centimeters', 'PaperPosition', [0 0 29.7 21], 'PaperOrientation', 'portrait');
subplot(1,2,1);
hold on;
for i=1:4
    hpatch = patch(issfs{i});
    set(hpatch,'FaceColor', colors{i}, 'EdgeColor', 'none');
    set(hpatch,'AmbientStrength', .2, 'SpecularStrength',.6, 'DiffuseStrength',.6);
end
plot3([5 256], [256 256], [zPlane zPlane], 'c', 'LineWidth', 2);
plot3([5 256], [5 5], [zPlane zPlane], 'c', 'LineWidth', 2);
plot3([5 5], [5 256], [zPlane zPlane], 'c', 'LineWidth', 2);
plot3([256 256], [5 256], [zPlane zPlane], 'c', 'LineWidth', 2);
% quiver3(150,30,30,1000/12,0,0, 'LineWidth', 2);
% quiver3(150,30,30,0,1000/12,0, 'LineWidth', 2);
% quiver3(150,30,30,0,0,1000/25, 'LineWidth', 2);
[x, y, z] = meshgrid(1, [1 256], [1 256]);
x = squeeze(x); y = squeeze(y); z = squeeze(z); 
surf(x, y, z, squeeze(raw(1,1:256,:)), 'Facecolor', 'texture', 'EdgeColor', 'none');
[x, y, z] = meshgrid([1 256], 1, [1 256]);
x = squeeze(x); y = squeeze(y); z = squeeze(z); 
surf(x, y, z, squeeze(raw(1:256,1,:)), 'Facecolor', 'texture', 'EdgeColor', 'none');
[x, y, z] = meshgrid([1 256], [1 256], 1);
x = squeeze(x); y = squeeze(y); z = squeeze(z); 
surf(x, y, z, squeeze(raw(1:256,1:256,1)), 'Facecolor', 'texture', 'EdgeColor', 'none');
caxis([50 200]);
colormap('gray');
lighting gouraud;
camproj perspective;
view(113,14);
hlight = camlight('headlight');
daspect([25 25 12]);
axis off;

subplot(1,2,2);
hold on;
imagesc(raw(256:-1:5,5:256,zPlane));
a = label2rgb(permute(trace(1:256, 1:256, zPlane), [2 1 3]), [1 0 0; 0 1 0; 0 0 1; 1 1 0], [1 1 1]);
h1 = imagesc(a(256:-1:5,5:256,:));
set(h1, 'AlphaData', .2);
axis equal;
axis off;

set(gcf,'PaperSize',fliplr(get(gcf,'PaperSize')));
print(gcf, 'C:\Users\mberning\Desktop\figure2.png', '-dpng', '-r900');

%% ADDED 2.12.2013: Mask and x-Affinity visualization
cnet.isoBorder = 2; [target, mask] = xyzMaskIso(cnet, trace);
%%
figure;
img = permute(target{2}(1:256,1:256,zPlane).*mask{2}(1:256,1:256,zPlane), [2 1 3]);
imagesc(img(5:256,5:256,:));
colormap('gray');
caxis([-1 1.5]);
axis equal;
axis off;

%% target or mask only

figure;
img = permute(mask{2}(1:256,1:256,zPlane), [2 1 3]);
imagesc(img(5:256,5:256,:));
colormap('gray');
caxis([-1 1.5]);
axis equal;
axis off;

