%% Intracellular vs extracellular hist cortex
clear all; clc;
% Intracellular stain
files = dir('D:\sync\trainingData\2013*.mat');
load(['D:\sync\trainingData\' files(1).name]);
a = raw(stack == 0); b = raw(stack ~= 0);
% Extracellular stain
load('C:\Users\mberning\Desktop\e2006_dense_final.mat');
raw = readKnossosRoi('E:\e2006\e2006_mag1\', '080823_e2006_mouseHRP_mag1', [5800 6000;3350 3550;1100 1300]);
[stack, bbox] = convertKleeTracingToLocalStack( KLEE_savedTracing );
stack = stack(51:150,51:150,51:150);
raw = raw(51:150,51:150,51:150);
c = raw(stack == 0); d = raw(stack ~= 0);

%%
x = 77.5:5:177.5;
c1 = hist(single(c),x);
d1 = hist(single(d),x);
a1 = hist(single(a),x);
b1 = hist(single(b),x);

figure;
subplot(2,2,1);
stairs(x, [c1; d1]', 'LineWidth', 2);
xlim([75 180]);
ylim([0 15e4]);
xlabel('voxel value');
ylabel('count');
subplot(2,2,2);
stairs(x, [a1; b1]', 'LineWidth', 2);
xlim([75 180]);
ylim([0 15e4]);
xlabel('voxel value');
ylabel('count');
set(findall(gcf,'-property','FontSize'),'FontSize',14)

