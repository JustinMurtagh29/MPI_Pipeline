both = load('D:\sync\activeTraining\classifierNewFirstTestBoth.mat');
skel = load('D:\sync\activeTraining\classifierNewFirstTestSkel.mat');

%% ARD results
figure;
plot(exp(skel.hyp.cov), 'r');
hold on;
plot(exp(both.hyp.cov), 'g');
set(gca, 'YScale', 'log');
plot([0 250], [10 10], 'b');
legend({'Skeletons' 'Both', 'Threshold proposal'});

%% Probability results
prob = sort(exp(skel.lp));
figure;
plot(prob, 'k', 'LineWidth', 2);
hold on;
idx1 = find(prob > .15, 1, 'first');
idx2 = find(prob > .95, 1, 'first');
p = patch([0 idx1 idx1 0],[0 0 1 1], 'r', 'EdgeColor', 'none');
set(p,'FaceAlpha',0.5);
p = patch([idx1 idx2 idx2 idx1],[0 0 1 1], 'b', 'EdgeColor', 'none');
set(p,'FaceAlpha',0.5);
p = patch([idx2 length(prob) length(prob) idx2],[0 0 1 1], 'g', 'EdgeColor', 'none');
set(p,'FaceAlpha',0.5);
legend({'Probability' 'Reject', 'Query', 'Accept'});
title('Skeleton based probabilities');

prob = sort(exp(both.lp));
figure;
plot(prob, 'k', 'LineWidth', 2);
hold on;
idx1 = find(prob > .15, 1, 'first');
idx2 = find(prob > .95, 1, 'first');
p = patch([0 idx1 idx1 0],[0 0 1 1], 'r', 'EdgeColor', 'none');
set(p,'FaceAlpha',0.5);
p = patch([idx1 idx2 idx2 idx1],[0 0 1 1], 'b', 'EdgeColor', 'none');
set(p,'FaceAlpha',0.5);
p = patch([idx2 length(prob) length(prob) idx2],[0 0 1 1], 'g', 'EdgeColor', 'none');
set(p,'FaceAlpha',0.5);
legend({'Probability' 'Reject', 'Query', 'Accept'});
title('Volume & Skeleton based probabilities');

%% Some numbers
label = generateWeightLabels();
threshold = 10;
clc;

display('Skeletons only training:');
display(['Rejected: ' num2str(sum(exp(skel.lp) < .15))]);
display(['Query: ' num2str(sum(exp(skel.lp) > .15 & exp(skel.lp) < .95))]);
display(['Accepted: ' num2str(sum(exp(skel.lp) > .95))]);
idx = exp(skel.hyp.cov) < threshold;
idx = find(idx);
display(['Relevant features:' num2str(length(idx))]);
for i=1:length(idx)
    display(label{idx(i)});
end

display('Both types training:');
display(['Rejected: ' num2str(sum(exp(both.lp) < .15))]);
display(['Query: ' num2str(sum(exp(both.lp) > .15 & exp(both.lp) < .95))]);
display(['Accepted: ' num2str(sum(exp(both.lp) > .95))]);
idx = exp(both.hyp.cov) < threshold;
idx = find(idx);
display(['Relevant features:' num2str(length(idx))]);
for i=1:length(idx)
    display(label{idx(i)});
end

%% Create histogram for each feature
label = generateWeightLabels();

for i=1:215
    visualizeEdgeFeaturesTrainingOnly( gtSkel, gtVol, label, i);
end
