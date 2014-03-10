%% Figure 4
load([param.dataFolder param.affSubfolder param.affMaps(map).name '.mat'], 'raw');
figure;
imagesc(raw(:,:,100));
colormap('gray');
axis off;
axis equal;
caxis([50 205]);
saveas(gcf, '/p/raw.pdf');
maximaC = [1 1 1 .5 .5];
for map=1:5
    load([param.dataFolder param.affSubfolder param.affMaps(map).name '.mat'], 'affX');
    figure;
    imagesc(affX(:,:,100));
    colormap('gray');
    axis off;
    axis equal;
    caxis([min(min(affX(:,:,100))) maximaC(map)]);
    saveas(gcf, ['/p/map' num2str(map) '.pdf']);
end

%% Figure 5
load([param.dataFolder param.outputSubfolder param.affMaps(map).name '/MorphRecon' num2str(r) '.mat']);
figure;
imagesc(v{1}(:,:,100))
colormap('gray');
axis off;
axis equal;
caxis([-0.1 1]);
saveas(gcf, '/p/morph.pdf');
seg = watershedSeg_v2_paper(v, .3, 150);
seg = seg{1};

%%
figure;
imagesc(fgm1(:,:,100))
colormap('gray');
axis off;
axis equal;
saveas(gcf, '/p/afterThres.pdf');

%%
figure;
imagesc(fgm1(:,:,100))
colormap('gray');
axis off;
axis equal;
saveas(gcf, '/p/afterVolume.pdf');

%%
figure;
load('segmentation/autoKLEE_colormap.mat');
autoKLEE_colormap = repmat(autoKLEE_colormap, 100, 1);
temp = label2rgb(segmentation{1}(:,:,100), autoKLEE_colormap, 'k');
imagesc(temp);
axis off;
axis equal;
saveas(gcf, '/p/seg.pdf');

%% Figure 6

load(param.cmSource);
autoKLEE_colormap = repmat(autoKLEE_colormap, 100, 1);
views = {2, 3, [90,0]};
for i=1:max(seg(:))
    obj = seg == i;
    issf{i} = isosurface(obj, .1);
end
display('Calculation done.');

%%
f = figure('Renderer', 'OpenGL', 'Visible', 'on' );
for currentView=1:length(views)
    subplot(1,3, currentView);
    hold on;
    for i=1:max(seg(:))
        k{i} = patch(issf{i});
        set(k{i}, 'FaceColor', autoKLEE_colormap(i,:), 'EdgeColor', 'none');
    end
    view(views{currentView});
    daspect([25 25 12]);
    grid on;
    alpha(.6);
    xlim([1 384]);
    ylim([1 384]);
    zlim([1 384]);
    camlight('headlight');
    lighting phong;
end

saveas(gcf, '/p/seg3D.tif');

