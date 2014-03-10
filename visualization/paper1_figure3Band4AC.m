clear all; clc
addpath(genpath('C:/code/CNN-Manuel/'));
addpath(genpath('C:/code/auxiliaryMethods/'));
data = findFiles('D:\cnetRetinaForPaper\CNNtrain\');

% Extract screened parameters
data2 = extractParameter(data);
data3 = data2(~cellfun(@isempty, {data2(:).name}));

% Sort according to date
dates = regexp({data3(:).name}, '-', 'split');
for i=1:length(dates)
    dates{i} = datenum([dates{i}{1:3}], 'ddmmmyyyy');
end
dates = [dates{:}];
[datesSorted, idxSort] = sort(dates, 'ascend');
data3 = data3(idxSort);

% Find training iterations of all CNN's
for i=1:length(data3)
	iterations(i) = length(data3(i).files).*data3(i).iterations;
end
clear i;

% Find 5 nets shown in paper
names = {data3(:).name};
idx(1) = find(strcmp(names, '14-May-2012-net957506'));
idx(2) = find(strcmp(names, '14-Jul-2012-net694915'));
idx(3) = find(strcmp(names, '14-Jul-2012-net777236'));
idx(4) = find(strcmp(names, '07-Sep-2012-net000133'));
idx(5) = find(strcmp(names, '07-Sep-2012-net618765'));

%% Generate supplementary parameter tables (OLD!!!)
% A = cat(1, data3(:).nrFM, data3(:).sizeFM, data3(:).nrHiddenLayer, data3(:).nrOutputs, data3(:).etaW, ...
%     data3(:).etaB, data3(:).tauW, data3(:).tauB, data3(:).batchSize, data3(:).iterations);
% xlswrite('C:\Users\mberning\Desktop\figureTemp\retinaParameterCNN.xls', reshape(A,532,10), 'B2:K533');
% xlswrite('C:\Users\mberning\Desktop\figureTemp\retinaParameterCNN.xls', {data3(:).name}', 'A2:A533');

%% Cortex Training Data Table
clear all; clc;
load D:\sync\trainingData\parameter.mat;
load D:\sync\trainingData\exclude.mat;
files = dir('D:\sync\trainingData\2013*.mat');
for i=1:length(files);
    if any(stacks(i).taskID == excludeTask)
        exclude(i) = 1;
    else
        exclude(i) = 0;
    end
end
for i=1:279
    load(['D:\sync\trainingData\' files(i).name]);
    voxelIntra(i) = sum(target(:) == 1);
    voxelExtra(i) = sum(target(:) == -1);
    voxelUnlabeled(i) = sum(target(:) == 0);
    meanValue(i) = mean(raw(:));
    raw = single(raw) - 122;
    stdValue(i) = std(raw(:));
end
for i=1:length(stacks)
    A{i,1} = stacks(i).taskID;
    A{i,2} = stacks(i).lastUpload;
    A{i,3} = stacks(i).tracer;
    A{i,4} = stacks(i).filename;
    A{i,5} = stacks(i).stackFile;
    A{i,6} = stacks(i).targetFile;
    A{i,7} = num2str(stacks(i).bbox(:)', '[%i,%i,%i;%i,%i,%i]');
    A{i,8} = num2str(stacks(i).bboxRaw(:)', '[%i,%i,%i;%i,%i,%i]');
    A{i,9} = exclude(i);
    A{i,10} = voxelIntra(i);
    A{i,11} = voxelExtra(i);
    A{i,12} = voxelUnlabeled(i);
    A{i,13} = meanValue(i);
    A{i,14} = stdValue(i);
end
xlswrite('C:\Users\mberning\Desktop\figureTemp\CNN_Cortex_Training_Data.xls', A, 'A2:N280');
xlswrite('C:\Users\mberning\Desktop\figureTemp\CNN_Cortex_Training_Data.xls', [fieldnames(stacks)' 'excluded' 'intracellular voxel' 'extracellular voxel' 'unlabelled voxel' 'mean value raw data' 'standard deviation raw data'], 'A1:N1');

%% Retina Training Data Table
[dataRaw, dataTrace] = getKleeStackList();
exclude = [6 26 62 99 122 123 124 149 153 171 174];
% for i=1:235
%     tracePresent = strcmp(['D:\sync\trainingData\retinaTransfer\tracing\e_k0563_ribbon_' num2str(i, '%.4i') '_fused.mat'], dataTrace);
%     if ~any(tracePresent);
%         fileToDelete = ['D:\sync\trainingData\retinaTransfer\raw\e_k0563_ribbon_' num2str(i, '%.4i') '_raw.mat'];
%         display(fileToDelete);
%         delete(fileToDelete);
%     end
% end
load([data3(idx(5)).directory data3(idx(5)).files{end}], 'cnet');
for i=1:length(dataTrace)
    bbox = zeros(3,2);
    load(dataRaw{i});
    load(dataTrace{i});
    [target, mask] = xyzMaskIso(cnet, kl_stack);
    voxelIntraTemp = 0;
    voxelExtraTemp = 0;
    voxelUnlabeledTemp = 0;
    for j=1:3
        temp = target{j}.*mask{j};
        voxelIntraTemp = voxelIntraTemp + sum(temp(:) == 1);
        voxelExtraTemp = voxelExtraTemp + sum(temp(:) == -1);
        voxelUnlabeledTemp = voxelUnlabeledTemp + sum(temp(:) == 0);
    end
    voxelIntra(i) = voxelIntraTemp;
    voxelExtra(i) = voxelExtraTemp;
    voxelUnlabeled(i) = voxelUnlabeledTemp;
    meanValue(i) = mean(kl_roi(:));
    kl_roi = single(kl_roi) - 122;
    stdValue(i) = std(kl_roi(:));
    A{i,1} = dataRaw{i};
    A{i,2} = dataTrace{i};
    if any(i == exclude)
        A{i,3} = 1;
    else
        A{i,3} = 0;
    end
    A{i,4} = voxelIntra(i);
    A{i,5} = voxelExtra(i);
    A{i,6} = voxelUnlabeled(i);
    A{i,7} = meanValue(i);
    A{i,8} = stdValue(i);
    temp2 = xlsread('I:\CortexConnectomics\Manuel\backup\20140129BigLinux\data\e_k0563\Matlab_Scripts\Fabian\Ribbon0001-0232_Knossoscoords_corrected.xls', ...
        'Tabelle1', ['B' num2str(str2double(dataRaw{i}(56:59))+1) ':' 'D' num2str(str2double(dataRaw{i}(56:59))+1)])';
    if ~isempty(temp2)
        bbox(:,1) = temp2;
        bbox(:,2) = bbox(:,1) + size(kl_roi)';
    end
    A{i,9} = num2str(bbox(:)', '[%i,%i,%i;%i,%i,%i]');
end
xlswrite('C:\Users\mberning\Desktop\figureTemp\CNN_Retina_Training_Data.xls', A, 'A2:I217');
xlswrite('C:\Users\mberning\Desktop\figureTemp\CNN_Retina_Training_Data.xls', {'dataRaw' 'dataTrace' 'excluded' 'intracellular voxel' 'extracellular voxel' 'unlabelled voxel' 'mean value raw data' 'standard deviation raw data' 'bounding box'}, 'A1:I1');

%% Supplementary table CNN Parameter cortex (TODO)
clear all; clc;
addpath(genpath('C:/code/CNN-cortex/'));
addpath(genpath('C:/code/auxiliaryMethods/'));
data = findFilesCortex('D:\cnetCortexForPaper\parameterSearch\');
data(cellfun(@isempty, data)) = [];
warning off; % Calling struct on a class seems to be warned about, needed here(?)

nrNets = 1;
for i=1:length(data)
    load(data{i});
    for iter=1:size(hyper.cnet,1)-1
        winnerIdx = [hyper.results{iter}.win(:).idx];
        for gpu=1:size(hyper.cnet,2)
            A{nrNets,1} = hyper.start;
            A{nrNets,2} = iter;
            A{nrNets,3} = gpu;
            % Rip out and strip apart classes xD
            run = struct(hyper.cnet(iter,gpu).run);
            run = rmfield(run, 'iterations');
            cnet = struct(hyper.cnet(iter,gpu));
            cnet = rmfield(cnet, 'run');
            cnet = rmfield(cnet, 'layer');
            cnetFields = fieldnames(cnet);
            runFields = fieldnames(run);
            % Write all parameter of one cnn instance in one line of A
            for j=1:length(cnetFields)
                if isa(cnet.(cnetFields{j}), 'function_handle')
                    A{nrNets,j+3} = char(cnet.(cnetFields{j}));
                elseif all(size(cnet.(cnetFields{j})) == [1 3])
                    A{nrNets,j+3} = num2str(cnet.(cnetFields{j}), '[%i,%i,%i]');
                elseif all(size(cnet.(cnetFields{j})) == [1 4])
                    A{nrNets,j+3} = num2str(cnet.(cnetFields{j}), '[%i,%i,%i,%i]');
                elseif all(size(cnet.(cnetFields{j})) == [1 5])
                    A{nrNets,j+3} = num2str(cnet.(cnetFields{j}), '[%i,%i,%i,%i,%i]');
                else
                    A{nrNets,j+3} = cnet.(cnetFields{j});
                end
            end
            for j=1:length(runFields)
                if isa(run.(runFields{j}), 'function_handle')
                    A{nrNets,j+length(cnetFields)+3} = char(run.(runFields{j}));
                elseif all(size(run.(runFields{j})) == [1 3])
                    A{nrNets,j+length(cnetFields)+3} = num2str(run.(runFields{j}), '[%i,%i,%i]');
                elseif all(size(run.(runFields{j})) == [1 4])
                    A{nrNets,j+length(cnetFields)+3} = num2str(run.(runFields{j}), '[%i,%i,%i,%i]');
                elseif all(size(run.(runFields{j})) == [1 5])
                    A{nrNets,j+length(cnetFields)+3} = num2str(run.(runFields{j}), '[%3.2e,%3.2e,%3.2e,%3.2e,%3.2e]');
                else
                    A{nrNets,j+length(cnetFields)+3} = run.(runFields{j});
                end
            end
            % This if is necessary if last net(s) of array did not finish 200 batches
            if gpu <= length(hyper.results{iter}.data) && ~isempty(hyper.results{iter}.data(gpu).err)
                A{nrNets,23} = mean(hyper.results{iter}.data(gpu).err(max(1,end-199):end)/1e6);
            else
                A{nrNets,23} = 'Inf';
            end
            A{nrNets,24} = any(gpu == winnerIdx);
            % Increase line counter for A
            nrNets = nrNets + 1;
        end
    end
end
xlswrite('C:\Users\mberning\Desktop\figureTemp\CNN_Cortex_Parameter.xls', A, 'A2:X3981');
xlswrite('C:\Users\mberning\Desktop\figureTemp\CNN_Cortex_Parameter.xls', ['startTime' 'iteration' 'GPU' cnetFields' runFields' 'mean squared voxel error last 200 stacks' 'keep/vary for next iteration'], 'A1:X1');

%% Supplementary Table Retina CNN parameter
clear par;
for i=1:length(data3)
    load([data3(i).directory data3(i).files{end}], 'cnet');
    run = struct(cnet.run);
    cnet = struct(cnet);
    cnet = rmfield(cnet, 'run');
    cnet = rmfield(cnet, 'layer');
    cnetFields = fieldnames(cnet);
    runFields = fieldnames(run);
    par.names = [cnetFields; runFields];
    for j=1:length(cnetFields)
        if isa(cnet.(cnetFields{j}), 'function_handle')
            par.A{i,j} = char(cnet.(cnetFields{j}));
        elseif all(size(cnet.(cnetFields{j})) == [1 3])
            par.A{i,j} = num2str(cnet.(cnetFields{j}), '[%i,%i,%i]');
        else
            par.A{i,j} = cnet.(cnetFields{j});
        end
    end
    for j=1:length(runFields)
        if isa(run.(runFields{j}), 'function_handle')
            par.A{i,j+length(cnetFields)} = char(run.(runFields{j}));
        elseif all(size(run.(runFields{j})) == [1 3])
            par.A{i,j+length(cnetFields)} = num2str(run.(runFields{j}), '[%i,%i,%i]');
        else
            par.A{i,j+length(cnetFields)} = run.(runFields{j});
        end
    end
end
xlswrite('C:\Users\mberning\Desktop\figureTemp\CNN_Retina_Parameter.xls', par.A, 'B2:AK533');
xlswrite('C:\Users\mberning\Desktop\figureTemp\CNN_Retina_Parameter.xls', {data3(:).name}', 'A2:A533');
xlswrite('C:\Users\mberning\Desktop\figureTemp\CNN_Retina_Parameter.xls', par.names', 'B1:AK1');

%% Figure 3B
values1 = [data3(:).nrFM];
values2 = [data3(:).nrHiddenLayer];
figure('Position', [1 41 1600 1084]);
plot(values1,values2, 'xr', 'MarkerSize', 20, 'LineWidth', 1.5);
xlim([0 35]);
ylim([0 10]);
xlabel('Feature maps');
ylabel('Hidden layer');
set(gca, 'box', 'off');
set(findall(gcf,'-property','LineWidth'),'LineWidth',1.5);
set(findall(gcf,'-property','FontSize'),'FontSize',28);
set(findall(gcf,'-property','Font'),'Font','Arial');
saveas(gcf, 'C:\Users\mberning\Desktop\figureTemp\fig3b.ai');

%% Figure 3C
values1 = [data3(:).sizeFM];
values2 = [data3(:).batchSize];
figure('Position', [1 41 1600 1084]);
plot(values1,values2, 'xr', 'MarkerSize', 20, 'LineWidth', 1.5);
set(gca, 'YScale', 'log');
set(gca, 'XScale', 'log');
xlabel('Filter size');
ylabel('Batch size');
set(gca, 'box', 'off');
set(findall(gcf,'-property','LineWidth'),'LineWidth',1.5);
set(findall(gcf,'-property','FontSize'),'FontSize',28);
set(findall(gcf,'-property','Font'),'Font','Arial');
saveas(gcf, 'C:\Users\mberning\Desktop\figureTemp\fig3c.ai');

%% Figure 3D
values1 = [data3(:).etaW];
values2 = [data3(:).etaB];
figure('Position', [1 41 1600 1084]);
plot(values1,values2, 'xr', 'MarkerSize', 20, 'LineWidth', 1.5);
set(gca, 'YScale', 'log');
set(gca, 'XScale', 'log');
xlabel('eta_w');
ylabel('eta_b');
set(gca, 'box', 'off');
set(findall(gcf,'-property','LineWidth'),'LineWidth',1.5);
set(findall(gcf,'-property','FontSize'),'FontSize',28);
set(findall(gcf,'-property','Font'),'Font','Arial');
saveas(gcf, 'C:\Users\mberning\Desktop\figureTemp\fig3d.ai');

%% Figure 3E
values1 = [data3(:).tauW];
values2 = [data3(:).tauB];
figure('Position', [1 41 1600 1084]);
plot(values1,values2, 'xr', 'MarkerSize', 20, 'LineWidth', 1.5);
set(gca, 'YScale', 'log');
set(gca, 'XScale', 'log');
xlabel('tau_w');
ylabel('tau_b');
set(gca, 'box', 'off');
set(findall(gcf,'-property','LineWidth'),'LineWidth',1.5);
set(findall(gcf,'-property','FontSize'),'FontSize',28);
set(findall(gcf,'-property','Font'),'Font','Arial');
saveas(gcf, 'C:\Users\mberning\Desktop\figureTemp\fig3e.ai');


%% Plot stuff for supp figure: application retina classifier to cortex 

% Load data to use
cortexData = readKnossosRoi('I:\CortexConnectomics\shared\cortex\2012-09-28_ex145_07x2\mag1\', '2012-09-28_ex145_07x2_mag1', [1501 1884; 1501 1884; 1501 1884] + [-18 17; -18 17; -8 7]);
retinaData = readKnossosRoi('I:\CortexConnectomics\shared\k0563_mag1\', '100527_k0563_mag1', [1664 2047; 1792 2175; 2176 2559] + [-18 17; -18 17; -8 7]);

% Load CNN (idx(5) = 618765 --> net for paper)
load([data3(idx(5)).directory data3(idx(5)).files{end}], 'cnet');
cnet.run.actvtTyp = @single;

% Normalize & fwd pass both stacks
retinaData = retinaData(:,:,81:120);
cortexData = cortexData(:,:,81:120);
retinaClass = cnet.fwdPass3D(cnet.normalizeStack(single(retinaData)));
cortexClass = cnet.fwdPass3D(cnet.normalizeStack(single(cortexData)));


%% Plot
plane = 22;

figure;
subplot(2,2,1);
imagesc(retinaData(1+(cnet.randOfConvn(1)-1)/2:end-(cnet.randOfConvn(1)-1)/2,1+(cnet.randOfConvn(2)-1)/2:end-(cnet.randOfConvn(2)-1)/2,plane+(cnet.randOfConvn(3)-1)/2));
colormap('gray');
axis equal; axis off;
title('raw retina');

subplot(2,2,2);
imagesc(single(retinaClass{6,1}(:,:,plane)));
colormap('gray');
caxis([-0.5 0.5]);
axis equal; axis off;
title('X affinity retina');

subplot(2,2,3);
imagesc(cortexData(1+(cnet.randOfConvn(1)-1)/2:end-(cnet.randOfConvn(1)-1)/2,1+(cnet.randOfConvn(2)-1)/2:end-(cnet.randOfConvn(2)-1)/2,plane+(cnet.randOfConvn(3)-1)/2));
colormap('gray');
axis equal; axis off;
title('raw cortex');

subplot(2,2,4);
imagesc(single(cortexClass{6,1}(:,:,plane)));
colormap('gray');
caxis([-0.5 0.5]);
axis equal; axis off;
title('X affinity cortex');
set(findall(gcf,'-property','FontSize'),'FontSize',16);

%% Plot stuff for figure 3F (make sure cnet is from idx(5) = 618765)
close all;
retinaData = readKnossosRoi('I:\CortexConnectomics\shared\k0563_mag1\', '100527_k0563_mag1', [1664 2047; 1792 2175; 2176 2559] + [-18 17; -18 17; -8 7]);
retinaData = retinaData(51:200,51:200,71:130);
retinaClass = cnet.fwdPass3D(cnet.normalizeStack(single(retinaData)));

%%
plane = 12;
for i=1:4
    figure;
    img = retinaClass{1,1}(1+17:end-18,1+17:end-18,plane+4+i);
    imagesc(img);
    caxis([-1.7 1.7]);
    colormap('gray');
    axis equal; axis off;
    saveas(gcf, ['C:\Users\mberning\Desktop\figureTemp\fig3f_raw' num2str(i) '.png']);
end

for i=[1,7,9]
    minVal = min(cnet.layer{2}.W{1,i}(:));
    maxVal = max(cnet.layer{2}.W{1,i}(:));
    figure;
    subplot(2,2,1);
    imagesc(cnet.layer{2}.W{1,i}(:,:,1));
    axis equal; axis off;
    caxis([minVal maxVal]);
    subplot(2,2,2);
    imagesc(cnet.layer{2}.W{1,i}(:,:,2));
    axis equal; axis off;
    caxis([minVal maxVal]);
    subplot(2,2,3);
    imagesc(cnet.layer{2}.W{1,i}(:,:,3));
    axis equal; axis off;
    caxis([minVal maxVal]);
    subplot(2,2,4);
    imagesc(cnet.layer{2}.W{1,i}(:,:,4));
    axis equal; axis off;
    caxis([minVal maxVal]);
    saveas(gcf, ['C:\Users\mberning\Desktop\figureTemp\fig3f_filter' num2str(i) '.png']);
end

for i=[1,7,9]
    minVal = min(retinaClass{2,i}(:));
    maxVal = max(retinaClass{2,i}(:));
    figure;
    img = retinaClass{2,i}(1+13:end-14,1+13:end-14,plane+6);
    imagesc(img);
    hold on;
    plot([1+3.5 1+3.5],[1+3.5 size(retinaClass{6,1},1)-3.5], '-r', 'LineWidth', 3);
    plot([1+3.5 size(retinaClass{6,1},1)-3.5],[1+3.5 1+3.5], '-r', 'LineWidth', 3);
    plot([size(retinaClass{6,1},1)-3.5 size(retinaClass{6,1},1)-3.5],[1+3.5 size(retinaClass{6,1},1)-3.5], '-r', 'LineWidth', 3);
    plot([1+3.5 size(retinaClass{6,1},1)-3.5],[size(retinaClass{6,1},1)-3.5 size(retinaClass{6,1},1)-3.5], '-r', 'LineWidth', 3);
    axis equal; axis off;
    caxis([minVal maxVal]);
    colormap('gray');
    saveas(gcf, ['C:\Users\mberning\Desktop\figureTemp\fig3f_fm' num2str(i) '.png']);
end

for i=1:3
    figure;
    imagesc(retinaClass{6,i}(:,:,plane));
    hold on;
    plot([1+17.5 1+17.5],[1+17.5 size(retinaClass{6,i},1)-17.5], '-r', 'LineWidth', 3);
    plot([1+17.5 size(retinaClass{6,i},1)-17.5],[1+17.5 1+17.5], '-r', 'LineWidth', 3);
    plot([size(retinaClass{6,i},1)-17.5 size(retinaClass{6,i},1)-17.5],[1+17.5 size(retinaClass{6,i},1)-17.5], '-r', 'LineWidth', 3);
    plot([1+17.5 size(retinaClass{6,i},1)-17.5],[size(retinaClass{6,i},1)-17.5 size(retinaClass{6,i},1)-17.5], '-r', 'LineWidth', 3);
    axis equal; axis off;
    caxis([-1 1]);
    colormap('gray');
    saveas(gcf, ['C:\Users\mberning\Desktop\figureTemp\fig3f_aff' num2str(i) '.png']);
end

%% Plot for figure 4B
retinaData = readKnossosRoi('I:\CortexConnectomics\shared\k0563_mag1\', '100527_k0563_mag1', [1664 2047; 1792 2175; 2176 2559] + [-18 17; -18 17; -8 7]);
retinaData = retinaData(:,:,71:130);

plane = 30;
for i=4:5
    load([data3(idx(i)).directory data3(idx(i)).files{end}]);
    cnet.run.actvtTyp = @single;

    % Normalize & fwd pass both stacks
    retinaClass = cnet.fwdPass3D(cnet.normalizeStack(single(retinaData)));

    % Plot
    figure;
    imagesc((retinaClass{6,1}(:,:,plane)+retinaClass{6,2}(:,:,plane)+retinaClass{6,3}(:,:,plane))/3); colormap('gray');
    axis equal; axis off;
    caxis([-1 1]); % Setting in paper: [-1 1]
    saveas(gcf, ['C:\Users\mberning\Desktop\figureTemp\fig4b' num2str(i) '.png']);
end

%% Calculate test error
clear retinaClass;
[filesRaw, filesTrace] = getKleeStackList();
load D:\retinaTrainingData\e_k0563_ribbon_0026_fused.mat;
load D:\retinaTrainingData\raw\e_k0563_ribbon_0026_raw.mat;
kl_roi = kl_roi(66:200,66:200,51:165);
kl_stack = kl_stack(66:200,66:200,51:165);
for dim=1:3
    inputPatch{dim} = 1:size(kl_stack,dim);
    outputPatch{dim} = inputPatch{dim}(ceil(cnet.randOfConvn(dim)/2:end - cnet.randOfConvn(dim)/2));
end
kl_stack = kl_stack(outputPatch{:});
[target, mask] = cnet.xyzMaskIso(kl_stack);

% Calculate error (takes quite some time)
err = cell(5,1);
for i=1:length(idx)
    display(['Starting ' num2str(i) 'th net error calculation:']);
    for j=1:length(data3(idx(i)).files)
        temp = [];
        tic;
        display(['File: ' num2str(j, '%.4i') '/' num2str(length(data3(idx(i)).files), '%.4i')]);
        load([data3(idx(i)).directory data3(idx(i)).files{j}], 'cnet');
        cnet.run.actvtTyp = @single;
        kl_class = cnet.fwdPass3D(cnet.normalizeStack(single(kl_roi)));
        for dim=1:3
            temp1 = .5.*(target{dim}-kl_class{6,dim}).^2;
            temp1 = temp1(mask{dim});
            temp(end+1:end+numel(temp1)) = single(temp1);
        end
        err{i}(j,1) = mean(temp);
        err{i}(j,2) = std(temp);
        toc;
    end
end
%save('D:\paperTemp\figure4C.mat');
%load('D:\paperTemp\figure4C.mat');

% Calculate error for random initalization
load([data3(idx(1)).directory data3(idx(1)).files{end}]);
cnet = cnet.init;
cnet.run.actvtTyp = @single;
kl_class = cnet.fwdPass3D(cnet.normalizeStack(single(kl_roi)));
temp = [];
for dim=1:3
    temp1 = .5.*(target{dim}-kl_class{6,dim}).^2;
    temp1 = temp1(mask{dim});
    temp(end+1:end+numel(temp1)) = single(temp1);
end
startError(1,1) = mean(temp);
startError(1,2) = std(temp);

% Load all old CNNs
clear cnet;
for i=1:5
    cnet{i} = load([data3(idx(i)).directory data3(idx(i)).files{end}]);
end

%% Figure 4C: CNN squared voxel error plot 
addpath(genpath('C:/code/CNN-Manuel/'));
load('D:\paperTemp\errorPlot.mat');
x{1} = (1:length(err{1}))*cnet{1}.cnet.run.maxIterMini;
x{2} = x{1}(end) + (1:length(err{2}))*cnet{2}.cnet.run.maxIterMini;
x{3} = x{1}(end) + (1:length(err{3}))*cnet{3}.cnet.run.maxIterMini;
x{4} = x{2}(end) + (1:length(err{4}))*cnet{4}.cnet.run.maxIterMini;
x{5} = x{2}(end) + (1:length(err{5}))*cnet{5}.cnet.run.maxIterMini;
names = {'CNN(1,1)','CNN(2,1)','CNN(2,2)','CNN(3,1)','CNN(3,2)'};
linestyles = {'xc','xg','xm','xy','xb'};

% Plot raw data
figure('Position', [1 41 1600 1084]);
for i=1:5
    subplot(5,3,(i-1)*3+1);
    errorbar(x{i}, err{i}(:,1), err{i}(:,2), linestyles{i});
    title(names{i});
    xlim([x{i}(1) x{i}(end)]);
    ylim([-1 4]);
    subplot(5,3,(i-1)*3+2);
    plot(x{i}, err{i}(:,1), linestyles{i});
    title([names{i} ' mean squared voxel error']);
    xlim([x{i}(1) x{i}(end)]);
    ylim([0 2]);
    subplot(5,3,(i-1)*3+3);
    plot(x{i}, err{i}(:,2), linestyles{i});
    title([names{i} ' std squared voxel error']);
    xlim([x{i}(1) x{i}(end)]);
    ylim([0 2]);
end
saveas(gcf, 'C:\Users\mberning\Desktop\figureTemp\fig4c1.ai');
saveas(gcf, 'C:\Users\mberning\Desktop\figureTemp\fig4c1.fig');

% Use sliding window (just on mean values)
lengthSlidingWindow = 101;
slidingWindow = @(x) conv(padarray(x, (lengthSlidingWindow - 1)/2, 'symmetric') , ones(lengthSlidingWindow,1)./lengthSlidingWindow, 'valid');
plot_variance = @(x,lower,upper,color) set(fill([x,x(end:-1:1)],[upper,fliplr(lower)],color),'EdgeColor',color,'FaceAlpha',0.2);
figure('Position', [1 41 1600 1084]);
plot([0 x{5}(end)], [startError(end,1) startError(end,1)], '-k');
hold on;
for i=1:5
    %plot_variance(x{i}, slidingWindow(err{i}(:,1))'+slidingWindow(err{i}(:,2))', slidingWindow(err{i}(:,1))'-slidingWindow(err{i}(:,2))', linestyles{i}(2));
    plot(x{i}, slidingWindow(err{i}(:,1)), linestyles{i});
end
xVert = (x{1}(end)+x{2}(1))/2;
plot([xVert xVert], [0 1], '--', 'Color', [.8 .8 .8]);
xVert = (x{2}(end)+x{4}(1))/2;
plot([xVert xVert], [0 1], '--', 'Color', [.8 .8 .8]);
legend([{'initalization'} names], 'Location', 'Best');
xlim([0 x{5}(end)]);
set(findall(gcf,'-property','MarkerSize'),'MarkerSize',5);
set(findall(gcf,'-property','LineWidth'),'LineWidth',1);
set(findall(gcf,'-property','FontSize'),'FontSize',14);
saveas(gcf, 'C:\Users\mberning\Desktop\figureTemp\fig4c2.eps');
saveas(gcf, 'C:\Users\mberning\Desktop\figureTemp\fig4c2.fig');

% Errorbar plot of used CNN versions (mean and std of mean within last #(slidingWindow) iterations of CNN)
figure('Position', [1 41 1600 1084]);
hold on;
for i=1:5
    errorbar(i, mean(err{i}(end-lengthSlidingWindow:end,1)), std(err{i}(end-lengthSlidingWindow:end,1)), linestyles{i});
end
set(findall(gcf,'-property','MarkerSize'),'MarkerSize',10);
set(findall(gcf,'-property','LineWidth'),'LineWidth',3);
set(findall(gcf,'-property','FontSize'),'FontSize',14);
set(gca, 'XTick', [1 2 3 4 5]);
set(gca, 'XTickLabel', names);
xlim([0 6]);
ylim([0 1.1]);
saveas(gcf, 'C:\Users\mberning\Desktop\figureTemp\fig4c3.ai');
saveas(gcf, 'C:\Users\mberning\Desktop\figureTemp\fig4c3.fig');

%% Plot Figure 5
addpath(genpath('C:/code/CNN-Manuel/'));
addpath(genpath('C:/code/auxiliaryMethods/'));
addpath(genpath('C:/code/segmentation/'));
class = retinaClass(6,1:3); clear retinaClass;
% mex -largeArrayDims -outdir segmentation/ segmentation/watershit_3D.cpp;
% Saved due to watershit crash, reproduced easily
%load('D:\paperTemp\figure5retina.mat'); % retina class from net 4
%completed as variable 'class' (only aff maps for memory purposes), not
%current version anymore!!!

% morphological reconstruction
rT = 1;
[x,y,z] = meshgrid(-rT:rT,-rT:rT,-rT:rT);
se = (x/rT).^2 + (y/rT).^2 + (z/rT).^2 <= 1;
for dir=1:3
    % Opening by reconstruction
    affEroded = imerode(class{dir}, se);
    affRecon = imreconstruct(affEroded, class{dir});
    % Closing by reconstruction
    affReconDilated = imdilate(affRecon, se);
    affReconRecon{dir} = imreconstruct(imcomplement(affReconDilated), imcomplement(affRecon));
end

% segmentation steps
affAll = (affReconRecon{1} + affReconRecon{2} + affReconRecon{3}) / 3;
fgm1 = affAll < .3;
fgm2 = bwareaopen(fgm1, 150, 26);
% load('D:\paperTemp\figure5retina2.mat'); % not up to date, still from
% using idx(4)
segmentation = watershed_3D(affReconRecon{1}, affReconRecon{2}, affReconRecon{3}, fgm2);

%%
plane = 102;

% raw plotten
figure;
imagesc(retinaData(1+(cnet.randOfConvn(1)-1)/2:end-(cnet.randOfConvn(1)-1)/2,1+(cnet.randOfConvn(2)-1)/2:end-(cnet.randOfConvn(2)-1)/2,plane+(cnet.randOfConvn(3)-1)/2));
colormap('gray');
axis off;
axis equal;
saveas(gcf, 'C:\Users\mberning\Desktop\figureTemp\fig5a.png');

% classification plotten
figure;
imagesc(class{3}(:,:,plane))
colormap('gray');
axis off;
axis equal;
caxis([-1 1]);
saveas(gcf, 'C:\Users\mberning\Desktop\figureTemp\fig5b.png');

% morphological reconstruction plotten
figure;
imagesc(affReconRecon{3}(:,:,plane));
colormap('gray');
axis off;
axis equal;
caxis([0 2]);
saveas(gcf, 'C:\Users\mberning\Desktop\figureTemp\fig5c.png');

% colormap to use
own_colormap = distinguishable_colors(max(segmentation(:)));
% old version
% load('segmentation/autoKLEE_colormap.mat');
% autoKLEE_colormap = repmat(autoKLEE_colormap, 100, 1);

% after threshold plot
figure;
temp = label2rgb(fgm1(:,:,plane) + fgm2(:,:,plane), [1 0 0; 1 1 1], 'k');
imagesc(temp);
axis off;
axis equal;
saveas(gcf, 'C:\Users\mberning\Desktop\figureTemp\fig5d.png');

% after removing small objects
figure;
temp = label2rgb(uint16(fgm2(:,:,plane)) .* segmentation(:,:,plane), own_colormap, 'k');
imagesc(temp);
axis off;
axis equal;
saveas(gcf, 'C:\Users\mberning\Desktop\figureTemp\fig5e.png');

% segmentation
figure;
temp = label2rgb(segmentation(:,:,plane), own_colormap, 'k');
imagesc(temp);
axis off;
axis equal;
saveas(gcf, 'C:\Users\mberning\Desktop\figureTemp\fig5f.png');

%% supplementary figure 1
load D:\retinaTrainingData\tracing\e_k0563_ribbon_0025_fused.mat;
load D:\retinaTrainingData\raw\e_k0563_ribbon_0025_raw.mat;
% Load CNN (idx(5) = 618765 --> net for paper)
load([data3(idx(5)).directory data3(idx(5)).files{end}], 'cnet');
cnet.run.actvtTyp = @single;
cnet.isoBorder = 1;
[target1, mask1] = cnet.xyzMaskIso(kl_stack);
cnet.isoBorder = 2;
[target2, mask2] = cnet.xyzMaskIso(kl_stack);
affMapsNames = {'X' 'Y' 'Z'};

%%
z = 150;
figure;
subplot(3,4,1);
imagesc(kl_roi(:,:,z));
axis equal; axis off;
title('raw data plane z');
subplot(3,4,5);
imagesc(kl_roi(:,:,z+1));
axis equal; axis off;
title('raw data plane z + 1');
colormap('gray');
map = distinguishable_colors(5);
subplot(3,4,2);
temp = label2rgb(kl_stack(:,:,z), map, 'k');
imagesc(temp);
axis equal; axis off;
title('annotation plane z');
subplot(3,4,6);
temp = label2rgb(kl_stack(:,:,z+1), map, 'k');
imagesc(temp);
axis equal; axis off;
title('annotation plane z + 1');
for i=1:3
    subplot(3,4,3+(i-1)*4);
    temp = target1{i}(:,:,z).*mask1{i}(:,:,z);
    temp(temp == -1) = 2;
    temp = label2rgb(temp, [1 1 1; 1 0 0], 'k');
    imagesc(temp);
    axis equal; axis off;
    title(['affinity ' affMapsNames{i} ' border = 1']);
end
for i=1:3
    subplot(3,4,4+(i-1)*4);
    temp = target2{i}(:,:,z).*mask2{i}(:,:,z);
    temp(temp == -1) = 2;
    temp = label2rgb(temp, [1 1 1; 1 0 0], 'k');
    imagesc(temp);
    axis equal; axis off;
    title(['affinity ' affMapsNames{i} ' border = 2']);
end
set(findall(gcf,'-property','FontSize'),'FontSize',16);