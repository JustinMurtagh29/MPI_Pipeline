function cnet = setParamsToActvtClass( cnet, varargin )
%SETPARAMSTOACTVTCLASS Set cnet params to actvtClass.
% INPUT varargin: Type/class conversion as an anonymous function, e.g.
%                 @single or @(x)gpuArray(single(x))
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

if nargin  > 1
    cnet.actvtClass = varargin{1};
end

for lyr = 2:cnet.layer
    cnet.W{lyr} = cnet.actvtClass(gather(cnet.W{lyr}));
    cnet.b{lyr} = cnet.actvtClass(gather(cnet.b{lyr}));
    if cnet.batchNorm(lyr)
        cnet.bn_beta{lyr} = cnet.actvtClass(gather(cnet.bn_beta{lyr}));
        cnet.bn_gamma{lyr} = cnet.actvtClass(gather(cnet.bn_gamma{lyr}));
        cnet.bn_muInf{lyr} = cnet.actvtClass(gather(cnet.bn_muInf{lyr}));
        cnet.bn_sig2Inf{lyr} = cnet.actvtClass(gather(cnet.bn_sig2Inf{lyr}));
    end
end

cnet.optimizer = cnet.optimizer.setParamsToActvtClass(cnet.actvtClass);

end

