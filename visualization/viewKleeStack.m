%% Load KLEE data
cd P:\code\not_in_repo\;
load P:\sync\toFermat\e_k0563_ribbon_0131_fused.mat; % Volume labeled data for ministack
load P:\sync\toFermat\e_k0563_ribbon_0131_raw.mat; % Raw data ministack
kl_roi = normalizeStack(double(kl_roi));
%posMinicube = [1571; 2486; 616]; % Position of KLEE minicube in KNOSSOS coord

%% Convert to stack?

kl_stack = convertKleeTracingToStack(KLEE_savedTracing, size(kl_roi), [0 0 0]);

%% Retrieve skeleton traces from disk 

load P:/data/skeletonLabel/skelTraces.mat;
% Unncomment in case one wants to retrieve skeletons from NML files insted
% of precalculated mat-File
% path{1} = 'P:\data\skeletonLabel\GCLek0563\';
% path{2} = 'P:\data\skeletonLabel\INLek0563\';
% skelTraces.GCL = readSkelTraceFolder(path{1});
% skelTraces.INL = readSkelTraceFolder(path{2});
% clear path;

%% Visualize minicube KLEE
load segmentation/autoKLEE_colormap.mat;
% colors = {'r' 'g' 'b' 'y' 'k' 'm' 'c'};
slices = {[] [] [1]};

%close all;
figure('Position', [-1599 1, 1600 1124], 'Renderer', 'OpenGL' );
hold on;
[k, issfs] = plotKLEEtracing(double(kl_stack), autoKLEE_colormap);
r = plotOriginalData(double(kl_roi), slices);
% plotSkelTracesMinicube(skelTraces.GCL, posMinicube, (size(kl_roi)-1)/2, 'c');
% plotSkelTracesMinicube(skelTraces.INL, posMinicube, (size(kl_roi)-1)/2, 'm');
view(37.5,20);
xlim([1 size(kl_roi,1)]);
ylim([1 size(kl_roi,1)]);
zlim([1 size(kl_roi,1)]);
set(gcf, 'Color', 'w');
daspect([25 25 12]);
axis off;
axis vis3d;
alpha(.7);
% zoom(1.5);
light = camlight('headlight');
lighting phong;

%% Export to Amira
exportSurfaceToAmira(issfs, 'amiraTest.mss');

%% Construct smaller cube
kl_stack = kl_stack(129:end-129, 129:end-129, 129:end-129);
kl_roi = kl_roi(129:end-129, 129:end-129, 129:end-129);

%% change KLEE tracing to fit redefined stack size

cellRem = [];
for i=1:size(KLEE_savedTracing.contours, 2)
    toRem = [];
    KLEE_savedTracing.contours{i} = KLEE_savedTracing.contours{i} - 128 * ones(size(KLEE_savedTracing.contours{i}));
    KLEE_savedTracing.contourList(i,2) = KLEE_savedTracing.contourList(i,2) - 128;
    for j=1:size(KLEE_savedTracing.contours{i}, 1)
        if sum(KLEE_savedTracing.contours{i}(j,:) <= 0) || sum(KLEE_savedTracing.contours{i}(j,:) > 256)
            toRem = [toRem j];
        end
    end
    KLEE_savedTracing.contours{i}(toRem,:) = [];
    if size(KLEE_savedTracing.contours{i},1) < 3
        cellRem = [cellRem i];
    end
end
KLEE_savedTracing.contours(cellRem') = [];
KLEE_savedTracing.contourList(cellRem',:) = [];
clear toRem cellRem i j;

%% Remove cell with ID 27 and 30

toRem = (KLEE_savedTracing.contourList(:,3) == 27);
KLEE_savedTracing.contours(toRem) = [];
KLEE_savedTracing.contourList(toRem,:) = [];
clear toRem;
toRem = (KLEE_savedTracing.contourList(:,3) == 30);
KLEE_savedTracing.contours(toRem) = [];
KLEE_savedTracing.contourList(toRem,:) = [];

%% Reomive objects resulting from overlap in contours

cellIDs = unique(KLEE_savedTracing.contourList(:,3));
for i=1:max(kl_stack(:))
    if i ~= cellIDs
        kl_stack(kl_stack == i) = 0; 
    end
end

%% find coordinates of each cell
cellIDs = unique(KLEE_savedTracing.contourList(:,3));
cellCoords = cell(max(KLEE_savedTracing.contourList(:,3)), 1);
for i=1:size(KLEE_savedTracing.contourList, 1)
    id = KLEE_savedTracing.contourList(i,3);
    cellCoords{id} = [cellCoords{id}; KLEE_savedTracing.contours{i}(2:end,:) - 128 * ones(size(KLEE_savedTracing.contours{i}(2:end,:)))];
end
cellCoords(cellfun('isempty', cellCoords)) = [];

struct.contourList = [ones(28,1),(1:28)',(1:28)'];
struct.contours = cellCoords;

%% 


%% Fly around
for i=1:36
   camorbit(10,0,'data')
   camlight(light, 'headlight')
   drawnow
end

%% Switch visibilty of cell indicated by second argument

switchVisibility(k, 2)

%% Camera Rotation

makeOriginalDataMovie(kl_roi);

%%

for i=1:30 
%    print(gcf, '-djpeg', '-r50', ['soonToBeMovie' num2str(i, '%04.0f') '.jpg']);
     campan(-0.5,0);
     drawnow
end
for i=1:30 
%    print(gcf, '-djpeg', '-r50', ['soonToBeMovie' num2str(i, '%04.0f') '.jpg']);
     campan(0.5,0);
     drawnow
end
