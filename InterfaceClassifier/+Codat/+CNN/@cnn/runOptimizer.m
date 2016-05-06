function cnet = runOptimizer( cnet, gradient )
%RUNOPTIMIZER Run the optimizer using the gadient to update cnet
%parameters.
% INPUT gradient: Gradient struct output from cnet.backprop.
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

paramsOld = cnet.param2Vec();
[paramsNew,optimizer] = cnet.optimizer.optimize(paramsOld,cnet.param2Vec(gradient));
cnet = cnet.vec2Param(paramsNew);
cnet.optimizer = optimizer;

end

