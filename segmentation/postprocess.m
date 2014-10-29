%% Load data
map = [4 5];
algo = [2 2];
r = [1 2];
par1 = [3 4];
par2 = [2 2];
load([param.dataFolder param.affSubfolder param.affMaps(map(1)).name '.mat'], 'raw');
for i=1:2
    load([param.dataFolder param.outputSubfolder param.affMaps(map(i)).name '/seg' num2str(r(i)) '-' num2str(algo(i)) '.mat']);
    seg{i} = v{par1(i),par2(i)};
end

%% Param struct reused
param.outputSubfolder = ['postprocess' '06-Nov-2012' filesep];
param.figureSubfolder = [param.outputSubfolder 'figures' filesep];
param.algo = {'v1'};
param.r = 0;
param.pR = {};
param.pR{1,1} = {1:2 500:500:3000};
param.affMaps = param.affMaps(1);

if ~exist([param.dataFolder param.outputSubfolder], 'dir')
    mkdir([param.dataFolder param.outputSubfolder]);
end

save([param.dataFolder param.outputSubfolder 'parameter.mat'], 'param');

%% Remove small objects and replot
cutoffs = 500:500:3000;
matlabpool 8;
segNew=cell(length(seg),length(cutoffs));
for j=1:length(seg)
    maxObjectId = max(seg{j}(:));
    sizeObj = hist(single(seg{j}(:)),single(maxObjectId+1));
    sizeObj(1) = [];
    if j==1
        load([param.dataFolder param.affSubfolder '09072012-618765.mat']);
    else
        load([param.dataFolder param.affSubfolder '09072012-000133.mat']);
    end
    parfor cut=1:length(cutoffs)
        temp = seg{j};
        for i=1:maxObjectId
            if sizeObj(i) < cutoffs(cut)
                temp(temp == i) = 0;
            end
        end
        segNew{j,cut} = watershed_3D_PP( affX, affY, affZ, temp);
    end
end
matlabpool close;
if ~exist([param.dataFolder param.outputSubfolder '05142012-957506/'], 'dir')
    mkdir([param.dataFolder param.outputSubfolder '05142012-957506/']);
end
parsave([param.dataFolder param.outputSubfolder '05142012-957506/seg1-1.mat'], segNew);

%% Alternative fill holes
segNew = seg;
matlabpool 2;
parfor j=1:length(seg)
    maxObjectId = max(seg{j}(:));
    for i=1:maxObjectId
        display([num2str(j) ': ' num2str(i, '%.3i')]);
        bw = segNew{j} == i;
        bw = imfill(bw, 'holes');
        segNew{j}(bw) = i;
    end
end
matlabpool close;
if ~exist([param.dataFolder param.outputSubfolder '05142012-957506/'], 'dir')
    mkdir([param.dataFolder param.outputSubfolder '05142012-957506/']);
end
parsave([param.dataFolder param.outputSubfolder '05142012-957506/seg1-1.mat'], segNew);

%% Evaluation
param.pR{1,1} = {1:6 1:6};
% Analysis
param = evalParameterSegmentation( param );
% Overview of performance of different segmentations
visualizeOverview( param );

%% Einzelne anschauen
map = 1;
algo = 1;
r = 1;
par1 = 2;
par2 = 6;
visualizeSingle(param, map, algo, r, par1, par2);

%% 
addpath('/home/mberning/code/KLEE');
KLEE_v4('stack', raw, 'stack_2', segNew{2,6});
KLEE_v4('stack', raw, 'stack_2', seg);
