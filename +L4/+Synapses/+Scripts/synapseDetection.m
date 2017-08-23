%% synapse detection pipeline from SynEM predictions
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

%% load data

p = Gaba.getSegParameters('ex145_ROI2017');
p = L4.updateParamsToNewestFiles(p);
graph = Seg.IO.loadGraph(p, true, 'end');