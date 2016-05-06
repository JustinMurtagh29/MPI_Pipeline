function cnet = setOptimizer( cnet, optimizer )
%SETOPTIMIZER Set optimizer of cnet.
% INPUT optimizer: A Codat.Optimizer object.
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

cnet.optimizer = optimizer.init(cnet.numParams);


end

