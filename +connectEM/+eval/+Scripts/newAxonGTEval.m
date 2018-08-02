% script to run the evaluation on the new_axon_gt axons
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

info = Util.runInfo();


%% get the overlaps between gt skeletons and agglos

p = Gaba.getSegParameters('ex145_ROI2017');
[skels, segIds] = connectEM.eval.getNewAxonGT();
for i = 1:10
    skels{i}.verbose = false;
end
m = load(p.agglo.axonAggloFile);
agglos = m.axons(m.indBigAxons);
ov = connectEM.eval.getNewAxonGTAggloOverlap(segIds, agglos);

[~, axFile] = fileparts(p.agglo.axonAggloFile);
outFile = fullfile(p.agglo.saveFolder, 'eval', ...
    sprintf('axon_gt_eval_%s.mat', axFile));

if ~exist(outFile, 'file')
    save(outFile, 'info', 'ov', 'skels', 'agglos');
end


%% everything to nml

ovSkels = connectEM.eval.newAxonGTOverlapsToSkel(skels, agglos, ov);

% write to tracings
for i = 1:10
    ovSkels{i}.write(sprintf('Axon%02d_ov1.nml', i));
end

