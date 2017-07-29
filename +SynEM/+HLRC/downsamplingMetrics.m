%% script for image comparison using synapse spatial frequencies
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

%% look at features used for synapse classification
m = load(['E:\workspace\data\backed\classifier\SynEM\' ...
    'SynEMPaperClassifier.mat'], 'classifier');
imp = m.classifier.ens.predictorImportance();
[~,idx] = sort(imp, 'descend');
fm = SynEM.getFeatureMap('paper');
names = fm.createNameStrings();
selFeat = false(length(idx), 1);
selFeat(idx(1:30)) = true;
fm.setSelectedFeat(selFeat);

%% MI for consecutive images

m = load('E:\workspace\sync\stackD_v9.mat', 'stackD');
stackD_9 = m.stackD;
m = load('E:\workspace\sync\stackD_v4.mat', 'stackD');
stackD_4 = m.stackD;
m = load('Cube67_raw.mat');
raw = m.raw;

[m_9, m_ref] = SynEM.HLRC.consMI(stackD_9, raw);
m_4 = SynEM.HLRC.consMI(stackD_4);
m_ref_cons = m_ref(m_ref < 1);

plot(m_4);
hold on
plot(m_9)
plot(m_ref_cons);
legend('v4','v9','L4 ref')
xlabel('z slice');
ylabel('mutual information')
a = gca;
a.FontSize = 24;
a.FontName = 'Arial';

%% syn spatial freq calculation

m = load('E:\workspace\sync\stackD_v9.mat', 'stackD');
stackD_9 = m.stackD;
m = load('E:\workspace\sync\stackD_v4.mat', 'stackD');
stackD_4 = m.stackD;
% m = load('Cube67_raw.mat');
% raw = m.raw;
% m = load('E:\workspace\allParametersFileServer.mat', 'p');
% p = m.p;
% bbox = p.local(67).bboxSmall;
% bbox = bsxfun(@plus, bbox(:,1), [ones(3, 1), size(stackD_9)']);
% raw = Seg.IO.loadRaw(p, bbox);
m = load('raw_ref.mat');
raw = m.raw;

raw = zscore(single(raw));
stackD_9 = zscore(single(stackD_9));
stackD_4 = zscore(single(stackD_4));

% object size in pixels
siz = [30 200]./11.24;
[ac_9, ac_ref] = SynEM.HLRC.synSpatialFreq(stackD_9, raw, siz, [1, 1, 2.5], 'pixelbpf');
[ac_4, ~] = SynEM.HLRC.synSpatialFreq(stackD_4, raw, siz, [1, 1, 2.5], 'pixelbpf');
