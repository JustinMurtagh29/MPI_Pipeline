%% Load raw data in two zoom level
addpath('P:\code\auxiliaryMethods\');
bbox = [1000 1500; 1000 1500; 1000 1250];
root = 'E:\e_k0563\';
zoom1 = 'k0563_mag1\';
zoom2 = 'k0563_mag2\';
normal = readKnossosRoi([root zoom1], '100527_k0563_mag1', bbox);
zoomed = readKnossosRoi([root zoom2], '100527_k0563_mag2', bbox/2);

%% Fix z bias similar to brainflight
normalF = zeros(size(normal).*[1 1 2]);
for i=1:2:length(normal)
   normalF(:,:,i) = normal(:,:,ceil(i/2));
   normalF(:,:,i+1) = normal(:,:,ceil(i/2));
end
zoomedF = zeros(size(zoomed).*[1 1 2]);
for i=1:2:length(zoomed)
   zoomedF(:,:,i) = zoomed(:,:,ceil(i/2));
   zoomedF(:,:,i+1) = zoomed(:,:,ceil(i/2));
end

%% Compare zoom levels

xyView = squeeze(normalF(129:256,129:256,250));
xzView = squeeze(normalF(129:256,250,129:256));
yzView = squeeze(normalF(250,129:256,129:256));
xyViewZ = squeeze(zoomedF(65:128,65:128,125));
xzViewZ = squeeze(zoomedF(65:128,125,65:128));
yzViewZ = squeeze(zoomedF(125,65:128,65:128));

figure;
subplot(2,3,1);
imagesc(xyView);
colormap('gray');
axis off;
subplot(2,3,2);
imagesc(xzView);
colormap('gray');
axis off;
subplot(2,3,3);
imagesc(yzView);
colormap('gray');
axis off;
subplot(2,3,4);
imagesc(xyViewZ);
colormap('gray');
axis off;
subplot(2,3,5);
imagesc(xzViewZ);
colormap('gray');
axis off;
subplot(2,3,6);
imagesc(yzViewZ);
colormap('gray');
axis off;

%% Duplicate all zoom voxels

for i=1:size(zoomedF,1)
    for j=1:size(zoomedF,2)
        for k=1:size(zoomedF,3)
            composed(2*i-1:2*i,2*j-1:2*j,2*k-1:2*k) = zoomedF(i,j,k);
        end
    end
end
intp = interp3(zoomedF,1);
%% Put into one image and plot
frame = [48 80];
figure;
for nr=1:500
    example = intp(129:256,129:256,nr);
    example2 = normalF(129:256,129:256,nr);
    finished = example;
    finished(frame(1):frame(2),frame(1):frame(2)) = example2(frame(1):frame(2),frame(1):frame(2));
    imagesc(finished);
    colormap('gray');
    axis off;
    drawnow;
    pause(.1);
end
