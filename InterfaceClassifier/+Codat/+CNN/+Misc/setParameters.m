function [ cnet ] = setParameters( cnet, W, b )
%SETPARAMETERS Set parameters of cnet to specified parameters.

for l = 2:length(W)
    cnet.W{l}(:,:,:,1:min(size(W{l},4),cnet.featureMaps(l-1)),1:min(size(W{l},5),cnet.featureMaps(l))) = W{l}(:,:,:,1:min(size(W{l},4),cnet.featureMaps(l-1)),1:min(size(W{l},5),cnet.featureMaps(l)));
    cnet.b{l}(1:min(length(b{l}),cnet.featureMaps(l))) = b{l}(1:min(length(b{l}),cnet.featureMaps(l)));
end


end
