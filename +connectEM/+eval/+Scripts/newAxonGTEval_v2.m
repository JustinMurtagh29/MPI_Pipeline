% script to evaluate the axon agglomeration using the new GT axons
%
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

info = Util.runInfo();


%% load axon agglo

p = Gaba.getSegParameters('ex145_ROI2017');
p.agglo.axonAggloFile = fullfile(p.agglo.saveFolder, ...
    'axons_19_a_partiallySplit_v2.mat');
axons = L4.Axons.getLargeAxons(p, true, true);


%% get the overlaps between gt skeletons and agglos


[skels, segIds] = connectEM.eval.getNewAxonGT();
for i = 1:10
    skels{i}.verbose = false;
end


m = load(p.agglo.axonAggloFile);
ov = connectEM.eval.getNewAxonGTAggloOverlap(segIds, axons);

[~, axFile] = fileparts(p.agglo.axonAggloFile);
outFile = fullfile(p.agglo.saveFolder, 'eval', ...
    sprintf('axon_gt_eval_%s.mat', axFile));

if ~exist(outFile, 'file')
    Util.ssave(outFile, info, ov, skels);
end


%% everything to nml

ovSkels = connectEM.eval.newAxonGTOverlapsToSkel(skels, agglos, ov);

% write to tracings
for i = 1:10
    ovSkels1{i}.write(sprintf('Axon%02d_ov1.nml', i));
end

