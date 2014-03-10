%% Load old data
clear all;
load P:\data\skeletonLabel\skeletonTraces.mat; % All skeleton data extracted with 2nd Cell of this skript
load P:\data\volumeLabel\rawTrace\e_k0563_ribbon_0120_fused.mat; % Volume labeled data for ministack
load P:\data\volumeLabel\rawTrace\e_k0563_ribbon_0120_raw.mat; % Raw data ministack
posMinicube = [1571; 2486; 616]; %aus AuswetungAmira.xls manuell anpassen an minicube
sizeMinicube = (size(kl_roi,1)-1)/2;

% Reshape KLEE tracing for visualization X-Y perm ute for KLEE done here
traceMatrices = convertContoursToVolumeLabel(KLEE_savedTracing);


%% Retrieve skeleton traces from disk 

path{1} = 'P:\data\skeletonLabel\GCLek0563';
path{2} = 'P:\data\skeletonLabel\INLek0563';

skelTraces.GCL = readSkelTraceFolder(path{1});
skelTraces.INL = readSkelTraceFolder(path{2});

%% Visualize skeletons
plotWhat = 'GCL';

figure('Position', [1600 -210, 1920 1120]);
hold on;
plotSkelTraces(skelTraces.(plotWhat));
view(3);

%% Visualize minicube KLEE

colors = {'r' 'g' 'b' 'y'};
slices = {[] [] [10:40:250]};

close all;
figure('Position', [1600 -210, 1920 1120], 'Renderer', 'OpenGL' );
hold on;
k = plotKLEEtracing(traceMatrices, colors);
r = plotOriginalData(double(kl_roi), slices);
% plotSkelTracesMinicube(skelTraces.GCL, posMinicube, sizeMinicube, 'c');
% plotSkelTracesMinicube(skelTraces.INL, posMinicube, sizeMinicube, 'm');
view(37.5,20);
xlim([1 size(kl_roi,1)]);
ylim([1 size(kl_roi,1)]);
zlim([1 size(kl_roi,1)]);
xlabel('x');
ylabel('y');
zlabel('z');
daspect([25 25 12]);
axis vis3d;
alpha(.8);
camlight;
lighting gouraud;
% for i=1:36
%    camorbit(10,0,'data')
%    drawnow
% end
