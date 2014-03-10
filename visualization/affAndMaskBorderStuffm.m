%% General Setup
close all;
xy = 110:160;
z = 150;
colors = distinguishable_colors(max(kl_stack(:)), [0 0 0]);
colors = [0 0 0; colors; 1 1 1];
figure('Position', [1 41 1600 1084]);

%% Stack
subplot(3,4,1);
imagesc(kl_stack(xy,xy,z));
colormap(colors);
caxis([0 max(kl_stack(:))+1]);
title('Training volume annotation (previous plane)');
axis equal; axis off;
subplot(3,4,5);
imagesc(kl_stack(xy,xy,z+1));
colormap(colors);
caxis([0 max(kl_stack(:))+1]);
title('Training volume annotation');
axis equal; axis off;

%% Targets
cnet.isoBorder = 1; [target1, mask1] = xyzMaskIso(cnet, kl_stack);
subplot(3,4,2);
imagesc(target1{1}(xy,xy,z).*mask1{1}(xy,xy,z));
title('x-affinity with target border = 1');
axis equal; axis off;
subplot(3,4,6);
imagesc(target1{2}(xy,xy,z).*mask1{2}(xy,xy,z));
title('y-affinity with target border = 1');
axis equal; axis off;
subplot(3,4,10);
imagesc(target1{3}(xy,xy,z).*mask1{3}(xy,xy,z));
title('z-affinity with target border = 1');
axis equal; axis off;

%% Targets with bigger border
cnet.isoBorder = 1; [target2, mask2] = xyzMaskIso(cnet, kl_stack);
target2{1} = imerode(target2{1}, ones(3,3,3));
target2{2} = imerode(target2{2}, ones(3,3,3));
target2{3} = imerode(target2{3}, ones(3,3,3));
subplot(3,4,3);
imagesc(target2{1}(xy,xy,z).*mask2{1}(xy,xy,z));
title('x-affinity with target border = 2');
axis equal; axis off;
subplot(3,4,7);
imagesc(target2{2}(xy,xy,z).*mask2{2}(xy,xy,z));
title('y-affinity with target border = 2');
axis equal; axis off;
subplot(3,4,11);
imagesc(target2{3}(xy,xy,z).*mask2{3}(xy,xy,z));
title('z-affinity with target border = 2');
axis equal; axis off;
