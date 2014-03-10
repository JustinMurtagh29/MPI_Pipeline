clear all; clc
addpath(genpath('C:/code/CNN-Manuel/'));
addpath(genpath('C:/code/auxiliaryMethods/'));
data = findFiles('D:\cnetRetinaForPaper\CNNtrain\');

%% Extract screened parameters
data2 = extractParameter(data);

%% Load network to use
idx = find(strcmp(data2.name(:), '14-May-2012-net957506'));
files = dir([data(idx).directory '*.mat']);
load([data(idx).directory files(end).name]);

%% Load data to use
cortexData = readKnossosRoi('I:\CortexConnectomics\shared\cortex\2012-09-28_ex145_07x2\mag1\', '2012-09-28_ex145_07x2_mag1', [1500 2000; 1500 2000; 1500 1520]);
retinaData = readKnossosRoi('I:\CortexConnectomics\shared\k0563_mag1\', '100527_k0563_mag1', [1500 2000; 1500 2000; 1500 1520]);

%% Apply CNN to data
retinaClass = cnet.fwdPass3D(normalizeStack(cnet, single(retinaData)));
cortexClass = cnet.fwdPass3D(normalizeStack(cnet, single(cortexData)));

%% Plot some stuff
figure; 
subplot(2,2,1);
imagesc(retinaData(1+(cnet.randOfConvn(1)-1)/2:end-(cnet.randOfConvn(1)-1)/2,1+(cnet.randOfConvn(2)-1)/2:end-(cnet.randOfConvn(2)-1)/2,1+(cnet.randOfConvn(3)-1)/2)); colormap('gray');
axis equal; axis off;
subplot(2,2,2);
imagesc(cortexData(1+(cnet.randOfConvn(1)-1)/2:end-(cnet.randOfConvn(1)-1)/2,1+(cnet.randOfConvn(2)-1)/2:end-(cnet.randOfConvn(2)-1)/2,1+(cnet.randOfConvn(3)-1)/2)); colormap('gray');
axis equal; axis off;
subplot(2,2,3);
imagesc(single(-retinaClass{6,3}(:,:,1))); colormap('gray'); caxis([-1.5 1.5]);
axis equal; axis off;
subplot(2,2,4);
imagesc(single(-cortexClass{6,3}(:,:,1))); colormap('gray'); caxis([-1.5 1.5]);
axis equal; axis off;
