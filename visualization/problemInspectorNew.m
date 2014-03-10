%% set this and run
clear all; clc;
id = '20130918T184016-00001';
subfolder = 'kdbLinearCov';
%subfolder = '';

%% Load data & preprocess
% Data writeKnowledgeDB
m = load(['D:\sync\problemInspector\' subfolder '\kDbResults' id '.mat']);
% Sort out possibleEnds
m.vShort = m.v;
for vp=1:length(m.vShort)
    for i=1:length(m.vShort(vp).missions)
        % ... that are not within viewport
        firstFrames = cellfun(@length, {m.vShort(vp).missions(i).possibleEnds(:).firstFrame});
        lastFrames = cellfun(@length, {m.vShort(vp).missions(i).possibleEnds(:).lastFrame});
        keep = firstFrames & lastFrames;
        m.vShort(vp).missions(i).possibleEnds = m.vShort(vp).missions(i).possibleEnds(keep);
        % ... that are farther away from the problem than [frames] frames
        frames = [40 40];
        distance = frames * 11.28;
        firstFrames = [m.vShort(vp).missions(i).possibleEnds(:).firstFrame] > distance(1);
        lastFrames = [m.vShort(vp).missions(i).possibleEnds(:).lastFrame] < -distance(2);
        toDel = firstFrames | lastFrames;
        m.vShort(vp).missions(i).possibleEnds = m.vShort(vp).missions(i).possibleEnds(~toDel);
    end
end
colors = [1 0 0; 0 1 0; distinguishable_colors(max(cellfun(@length, {m.vShort(2).missions(:).possibleEnds}))-1, [1 0 0; 0 1 0])];

%% Data fromGraphToDB
load(['D:\sync\problemInspector\' subfolder '\kDbBefore' id '.mat']);
% Weight label for feature visualization (needs update if features are updated)
label = generateWeightLabels();
% Find wrong classification in training data
isConnected = zeros(size(edges,1),1);
for i=1:size(edges,1)
    % Compare only with skeleton training data
    edgeFound = all(bsxfun(@eq, train.edges(1:23491,:), edges(i,:)),2);
    if any(edgeFound)
        allLabels = train.y(edgeFound);
        isConnected(i) = unique(allLabels);
    end
end
logicalMerge = prob > .95;
logicalDiscard = prob < .15;
% Feedback evaluation (type 1 & 2 error etc)
lowerThres = 0:1:100;
upperThres = 0:1:100;
for i=1:length(lowerThres)
    pos = prob < (lowerThres(i)/100);
    neg = prob >= (lowerThres(i)/100);
    lTP(i) = sum(isConnected(pos) == -1)./length(isConnected)*100;
    lFP(i) = sum(isConnected(pos) == 1)./length(isConnected)*100;
    lUP(i) = sum(isConnected(pos) == 0)./length(isConnected)*100;
    lTN(i) = sum(isConnected(neg) == 1)./length(isConnected)*100;
    lFN(i) = sum(isConnected(neg) == -1)./length(isConnected)*100;
    lUN(i) = sum(isConnected(neg) == 0)./length(isConnected)*100;
    lPos(i) = sum(pos)./length(isConnected)*100; % count only labeled points in this case
end
for i=1:length(upperThres)
    pos = prob > (upperThres(i)/100);
    neg = prob <= (upperThres(i)/100);
    uTP(i) = sum(isConnected(pos) == 1)./length(isConnected)*100;
    uFP(i) = sum(isConnected(pos) == -1)./length(isConnected)*100;
    uUP(i) = sum(isConnected(pos) == 0)./length(isConnected)*100;
    uTN(i) = sum(isConnected(neg) == -1)./length(isConnected)*100;
    uFN(i) = sum(isConnected(neg) == 1)./length(isConnected)*100;
    uUN(i) = sum(isConnected(neg) == 0)./length(isConnected)*100;
    uPos(i) = sum(pos)./length(isConnected)*100;
end

%% Global Statistics on batch of missions
visualizeRenderingLength(m.vShort, id);

%% Mission Visualization, Statistics, Videos to debug levelcreator
for i=9:50
    visualizeTasksNewNew(m.seg, m.vShort, m.par, id, i, colors);
    %arbitraryReslice(m.raw, m.seg, m.mito, m.vShort, m.par, id, i, colors);
end

%% plot classification & threshold performance
close all;
[probTemp, idx] = sort(prob);
figure('Units', 'centimeters', 'Position', [0 0 29.7 21], ...
'PaperType', 'A4', 'PaperUnits', 'centimeters', 'PaperPosition', [0 0 29.7 21], 'PaperOrientation', 'portrait');
plot(probTemp, 'k', 'LineWidth', 2);
hold on;
idx1 = find(probTemp > .15, 1, 'first');
idx2 = find(probTemp > .95, 1, 'first');
p = patch([0 idx1 idx1 0],[0 0 1 1], 'r', 'EdgeColor', 'none');
set(p,'FaceAlpha',0.5);
p = patch([idx1 idx2 idx2 idx1],[0 0 1 1], 'b', 'EdgeColor', 'none');
set(p,'FaceAlpha',0.5);
p = patch([idx2 length(probTemp) length(probTemp) idx2],[0 0 1 1], 'g', 'EdgeColor', 'none');
set(p,'FaceAlpha',0.5);
x = 1:length(probTemp);
f(:,1) = exp(feval(likfunc,[],ones(size(c)),c(idx)-2*sqrt(d(idx)), d(idx)));
f(:,2) = exp(feval(likfunc,[],ones(size(c)),c(idx)+2*sqrt(d(idx)), d(idx)));
plot(f(:,1), 'c');
plot(f(:,2), 'c');
xlim([0 size(probTemp,1)]);
xlabel('Edges sorted according to probability');
ylabel('Probability acooriding to classifier');
legend({'Probability', 'Reject', 'Query', 'Accept', 'Confidence'}, 'Location', 'Northwest');
title('Classification statistics for this batch');
saveas(gcf, 'C:\Users\mberning\Desktop\classifier\plot1.fig');
set(gcf,'PaperSize',fliplr(get(gcf,'PaperSize')));
print(gcf, 'C:\Users\mberning\Desktop\classifier\plot1.pdf', '-dpdf');

%% Type 1 & 2 error visualization for both classifications
figure('Units', 'centimeters', 'Position', [0 0 29.7 21], ...
'PaperType', 'A4', 'PaperUnits', 'centimeters', 'PaperPosition', [0 0 29.7 21], 'PaperOrientation', 'portrait');
subplot(2,1,1);
plot([lTP; lFP; lUP; lTN; lFN; lUN]', 'LineWidth', 2);
set(gca, 'XTick', 1:5:length(lowerThres));
set(gca, 'XTickLabel', lowerThres(1:5:length(lowerThres)));
xlim([1 length(lowerThres)]);
legend({'true positive' 'false positive' 'unlabeled positive' 'true negative' 'false negative' 'unlabeled negative'});
title('lowerThreshold: Is edge not connected?');
xlabel('threshold');
ylabel('percentage of classified edges');
subplot(2,1,2);
plot([uTP; uFP; uUP; uTN; uFN; uUN]', 'LineWidth', 2);
set(gca, 'XTick', 1:5:length(upperThres));
set(gca, 'XTickLabel', upperThres(1:5:length(upperThres)));
xlim([1 length(upperThres)]);
legend({'true positive' 'false positive' 'unlabeled positive' 'true negative' 'false negative' 'unlabeled negative'});
title('upperThreshold: Is edge connected?');
xlabel('threshold');
ylabel('percentage of classified edges');
saveas(gcf, 'C:\Users\mberning\Desktop\classifier\plot2.fig');
set(gcf,'PaperSize',fliplr(get(gcf,'PaperSize')));
print(gcf, 'C:\Users\mberning\Desktop\classifier\plot2.pdf', '-dpdf');

%% Precision recall curve
colors = jet(256);
precision = lTP ./ (lTP + lFP);
precision(1) = 1;
recall = lTP ./ (lTP + lFN);
p = uTP ./ (uTP + uFP);
p(end) = 1;
r = uTP ./ (uTP + uFN);
figure('Units', 'centimeters', 'Position', [0 0 29.7 21], ...
'PaperType', 'A4', 'PaperUnits', 'centimeters', 'PaperPosition', [0 0 29.7 21], 'PaperOrientation', 'portrait');
hold on;
for i=1:length(recall)
    scatter(recall(i), precision(i), 50, lPos(i), 'x', 'LineWidth', 2);
    scatter(r(i), p(i), 50, uPos(i), 'o', 'LineWidth', 2); 
end
colormap(colors);
colorbar;
title('Color indicates percent of edges discarded');
xlim([0 1]);
ylim([0 1]);
caxis([0 max(lPos)]);
xlabel('recall');
ylabel('precision');
legend('lowerThreshold', 'upperThreshold', 'Location', 'SouthWest');
saveas(gcf, 'C:\Users\mberning\Desktop\classifier\plot2_2.fig');
set(gcf,'PaperSize',fliplr(get(gcf,'PaperSize')));
print(gcf, 'C:\Users\mberning\Desktop\classifier\plot2_2.pdf', '-dpdf');

%% Accuracy sensitivity curve
sensitivity = lTN ./ (lTN + lFP);
sensitivity(1) = 1;
accuracy = (lTP + lTN) ./ (lTP + lTN + lFP + lFN);
s = uTN ./ (uTN + uFP);
s(end) = 1;
a = (uTP + uTN) ./ (uTP + uTN + uFP + uFN);
figure('Units', 'centimeters', 'Position', [0 0 29.7 21], ...
'PaperType', 'A4', 'PaperUnits', 'centimeters', 'PaperPosition', [0 0 29.7 21], 'PaperOrientation', 'portrait');
hold on;
for i=1:length(sensitivity)
    scatter(sensitivity(i), accuracy(i), 50, lPos(i), 'x', 'LineWidth', 2);
    scatter(s(i), a(i), 50, uPos(i), 'o', 'LineWidth', 2);
end
colormap(colors)
colorbar;
title('Accuracy over Sensitivity (Color = % of edges discarded due to upper or lower threshold, 34k without any drops)');
xlim([0 1]);
ylim([0 1]);
caxis([0 max(lPos)]);
xlabel('sensitivity');
ylabel('accuracy');
legend('lowerThreshold', 'upperThreshold', 'Location', 'SouthWest');
saveas(gcf, 'C:\Users\mberning\Desktop\classifier\plot2_3.fig');
set(gcf,'PaperSize',fliplr(get(gcf,'PaperSize')));
print(gcf, 'C:\Users\mberning\Desktop\classifier\plot2_3.pdf', '-dpdf');

%% Special plot for Moritz Wrong feedback over queries
figure('Units', 'centimeters', 'Position', [0 0 29.7 21], ...
'PaperType', 'A4', 'PaperUnits', 'centimeters', 'PaperPosition', [0 0 29.7 21], 'PaperOrientation', 'portrait');
feedbackRate = (lFP./(lTP+lFP)*100)';
feedbackRate(1) = 0;
plot(100 - lPos, feedbackRate, 'LineWidth', 2);
title('lowerThreshold: Querry number vs. wrong feedback tradeoff');
xlabel('number queries remaining (positives over all edges) [%]');
ylabel('wrong feedback (true over all positives) [%]');
saveas(gcf, 'C:\Users\mberning\Desktop\classifier\plot3.fig');
set(gcf,'PaperSize',fliplr(get(gcf,'PaperSize')));
print(gcf, 'C:\Users\mberning\Desktop\classifier\plot3.pdf', '-dpdf');

%% Automatic Relevance Detection results
close all;
figure('Units', 'centimeters', 'Position', [0 0 29.7 21], ...
'PaperType', 'A4', 'PaperUnits', 'centimeters', 'PaperPosition', [0 0 29.7 21], 'PaperOrientation', 'portrait');
thres = 7;
plot(exp(hyp.cov), 'r');
hold on;
set(gca, 'YScale', 'log');
plot([0 215], [thres thres], 'b');
idx = find(exp(hyp.cov) < thres);
items = {'Hyperparameter ARD', 'Threshold proposal'};
for i=1:length(idx)
    plot([idx(i) idx(i)], [.1 1000], 'c');
    items{end+1} = label{idx(i)};
end
xlim([0 215]);
legend(items, 'Location', 'BestOutside');
saveas(gcf, 'C:\Users\mberning\Desktop\classifier\plot4.fig');
set(gcf,'PaperSize',fliplr(get(gcf,'PaperSize')));
print(gcf, 'C:\Users\mberning\Desktop\classifier\plot4.pdf', '-dpdf');

%% Create histogram for each feature
for i=145
    visualizeEdgeFeaturesNewest( train, test, label, i);
end

%% Visualize probability distributions for each 
figure('Units', 'centimeters', 'Position', [0 0 29.7 21], ...
'PaperType', 'A4', 'PaperUnits', 'centimeters', 'PaperPosition', [0 0 29.7 21], 'PaperOrientation', 'portrait');
subplot(2,1,1);
hist(prob(logicalMerge & logicalError),100);
title('Predicted probabilities: False Merger');
subplot(2,1,2);
hist(prob(logicalDiscard & logicalError),100);
title('Predicted probabilities: False Discards');
saveas(gcf, 'C:\Users\mberning\Desktop\classifier\plot5.fig');
set(gcf,'PaperSize',fliplr(get(gcf,'PaperSize')));
print(gcf, 'C:\Users\mberning\Desktop\classifier\plot5.pdf', '-dpdf');
%% Generate movies for each false merger/discard and according feature visualization
idx = find((logicalMerge & (isConnected == -1)) | logicalDiscard & (isConnected == 1));
for i=191:length(idx)
    figure('Units', 'centimeters', 'Position', [0 0 29.7 21], ...
    'PaperType', 'A4', 'PaperUnits', 'centimeters', 'PaperPosition', [0 0 29.7 21], 'PaperOrientation', 'portrait');
    visualizeEdgeFeaturesForErrors(test.X, label, idx(i), logicalMerge, logicalDiscard, isConnected);
    if logicalMerge(idx(i))
        prefix = 'merge';
    end
    if logicalDiscard(idx(i))
        prefix = 'discard';
    end
    saveas(gcf, ['C:\Users\mberning\Desktop\classifier\errors\' prefix num2str(idx(i)) '.fig']);
    set(gcf,'PaperSize',fliplr(get(gcf,'PaperSize')));
    print(gcf, ['C:\Users\mberning\Desktop\classifier\errors\' prefix num2str(idx(i)) '.pdf'], '-dpdf');
    close all;
    arbitraryResliceForErrors(raw, seg, edges(idx(i),:), colors, idx(i), prefix);
    display(num2str(i));
end

%% Plot a nice visualization for constructing the supervoxel graph
for i=9:60
    visualizeTasksSupervoxel(m.seg, m.vShort, m.par, i);
end

