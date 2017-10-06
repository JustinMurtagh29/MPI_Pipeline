% script to run the evaluation on the new_axon_gt axons
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

%% get the overlaps between gt skeletons and agglos

p = Gaba.getSegParameters('ex145_ROI2017');
[skels, segIds] = connectEM.eval.getNewAxonGT();
for i = 1:10
    skels{i}.verbose = false;
end
m = load(p.agglo.axonAggloFile);
agglos = m.axons(m.indBigAxons);
ov = connectEM.eval.getNewAxonGTAggloOverlap(segIds, agglos);

%% everything to nml

ovSkels1 = connectEM.eval.newAxonGTOverlapsToSkel(skels, agglos, ov);
ovSkels2 = connectEM.eval.newAxonGTOverlapsToSkel(skels, agglos, ov, 2);
