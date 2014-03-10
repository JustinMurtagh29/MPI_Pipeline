clear all;
clc;
addpath('visualization');
[dataRaw, dataTrace] = getKleeStackList();
load segmentation/autoKLEE_colormap.mat;
slices = {[1] [1] [1]};


for i=1:length(dataRaw)
    figure('Position', [-1599 1, 1600 1124], 'Renderer', 'OpenGL', 'Visible', 'off');
    load(dataRaw{i});
    load(dataTrace{i});
    hold on;
    plotKLEEtracing(double(kl_stack), autoKLEE_colormap);
    plotOriginalData(double(kl_roi), slices);
    xlim([1 size(kl_roi,1)]);
    ylim([1 size(kl_roi,1)]);
    zlim([1 size(kl_roi,1)]);    
    daspect([25 25 12]);
    axis off;
    view(110,25);
    camlight('headlight');
    lighting phong;
    set(gcf, 'Color', 'w');
    set(gcf,'PaperPositionMode', 'manual', 'PaperUnits','centimeters', ...
    'Paperposition',[1 1 28 20], 'PaperSize', [29.2 21])
    drawnow;
    print(gcf, '-dpdf', '-r300', ['C:\Users\mberning\Desktop\figuresPaper\smallRibbons\ribbon' num2str(i, '%.3i')]);
    close(gcf);
end

%% Combine TIFF (use Photoshop Image Processor)

files = dir('C:\Users\mberning\Desktop\figuresPaper\smallRibbons\TIFF\*.tif');
figure('Position', [-1599 1, 1600 1124], 'Renderer', 'OpenGL', 'Visible', 'off');
for i=1:length(files)
    subaxis(11, 20, i, 'Spacing', 0.001, 'Padding', 0, 'Margin', 0);
    a = imread(['C:\Users\mberning\Desktop\figuresPaper\smallRibbons\TIFF\' files(i).name]);
    imagesc(a);
    axis tight;
    axis off;
end
set(gcf, 'Color', 'w');
set(gcf,'PaperPositionMode', 'auto', 'PaperUnits','centimeters', ...
    'PaperSize', [50 50]);
drawnow;
print(gcf, '-dpdf', '-r600', 'C:\Users\mberning\Desktop\figuresPaper\smallRibbon');


%% Plot one ribbon for big picture
load(dataTrace{130});
load(dataRaw{130});
load segmentation/autoKLEE_colormap.mat;
% colors = {'r' 'g' 'b' 'y' 'k' 'm' 'c'};
slices = {[1] [1] [1]};

figure('Position', [-1599 1, 1600 1124], 'Renderer', 'OpenGL' );
hold on;
[k, issfs] = plotKLEEtracing(double(kl_stack), autoKLEE_colormap);
r = plotOriginalData(double(kl_roi), slices);
xlim([1 size(kl_roi,1)]);
ylim([1 size(kl_roi,1)]);
zlim([1 size(kl_roi,1)]);
set(gcf, 'Color', 'w');
daspect([25 25 12]);
axis off;
view(110,25)
light = camlight('headlight');
lighting phong;

set(gcf,'PaperPositionMode', 'manual', 'PaperUnits','centimeters', ...
    'Paperposition',[1 1 28 20], 'PaperSize', [29.2 21])
drawnow;
print(gcf, '-dpdf', '-r300', 'C:\Users\mberning\Desktop\figuresPaper\bigRibbon');

%% Same image stichting for Figure 7 as Figure 2

%% Combine TIFF (use Photoshop Image Processor)

files = dir('C:\Users\mberning\Desktop\New folder\TIFF\*.tif');
figure('Position', [-1599 1, 1600 1124], 'Renderer', 'OpenGL', 'Visible', 'off');
for i=1:length(files)
    subaxis(11, 20, i, 'Spacing', 0.001, 'Padding', 0, 'Margin', 0);
    a = imread(['C:\Users\mberning\Desktop\New folder\TIFF\' files(i).name]);
    imagesc(a);
    axis tight;
    axis off;
end
set(gcf, 'Color', 'w');
set(gcf,'PaperPositionMode', 'auto', 'PaperUnits','centimeters', ...
    'PaperSize', [50 50]);
drawnow;
print(gcf, '-dpdf', '-r600', 'C:\Users\mberning\Desktop\figuresPaper\smallObjChains');


