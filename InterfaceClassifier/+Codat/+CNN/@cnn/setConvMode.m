function [ cnet ] = setConvMode( cnet, mode )
%SETCONVMODE Set mode for convolutions and adapt filters.
% INPUT cnet: Codat.CNN object.
%       mode: String specifying convolution mode
%           'fast': Fastest mode, high memory requirements
%           'memory1': Medium speed, medium memory requirements
%           'memory2:' Slowest speed, lowest memory requirements (usually
%                      not necessary to use it).
% NOTE This function is necessary since the two memory modes need the
% filter to be explicitly strided while the fast mode uses fft properties
% to introduce the stride in fourier space.
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

switch mode
    case {'fft1','fft2'}
        if strcmp('convn',cnet.convAlg)
            for lyr = 2:cnet.layer
                cnet.W{lyr} = reshape(cnet.W{lyr}(repmat(cnet.Wmask{lyr},[1 1 1 cnet.featureMaps(lyr - 1) cnet.featureMaps(lyr)])),[cnet.filterSize{lyr} cnet.featureMaps(lyr - 1) cnet.featureMaps(lyr)]);
            end
        end
        cnet.convAlg = mode;
    case 'convn'
        if any(strcmp({'fft1','fft2'},cnet.convAlg))
            for lyr = 2:cnet.layer
                cnet.W{lyr} = cnet.sparseKernel(cnet.W{lyr},cnet.d(lyr - 1,:));
            end
        end
        cnet.convAlg = mode;
    otherwise
        error('Unknown mode.');
end

end

