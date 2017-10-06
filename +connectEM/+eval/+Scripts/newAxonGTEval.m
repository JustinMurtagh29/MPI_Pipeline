% script to run the evaluation on the new_axon_gt axons
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

%% get the overlaps between gt skeletons and agglos

p = Gaba.getSegParameters('ex145_ROI2017');
[skels, segIds] = connectEM.eval.getNewAxonGT();
m = load(p.agglo.axonAggloFile);
agglos = m.axons(m.indBigAxons);
ov = connectEM.eval.getNewAxonGTAggloOverlap(segIds, agglos);

%% everything to nml

ovSkels = connectEM.eval.newAxonGTOverlapsToSkel(skels, agglos, ov);