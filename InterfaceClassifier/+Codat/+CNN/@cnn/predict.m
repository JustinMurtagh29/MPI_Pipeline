function [ prediction ] = predict( cnet, input )
%PREDICT Prediction for input cube.
% The difference to forward pass is that all intermediate results for
% backpropagation are discarded to have more memory available.
% INPUT input: 4d cube where the first three dimensions are the data and
%              the forth dimension is the number of input channels.
% OUTPUT prediction: 4d cube of voxel labels.
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

%initialize cells
input = cnet.actvtClass(input);
activity = input;
activityShortcut = cell(cnet.layer,1);

for lyr = 2:cnet.layer
    
    %conv layer
    activity = cnet.convLayerFwd(activity, lyr);
    
    %mp layer
    if cnet.maxPool(lyr)
        activity = cnet.maxPooling(activity,lyr);
    end
    
    %shortcut in
    if cnet.shortcut(lyr) > 0
        activity = activity + cnet.cropActivation(activityShortcut{lyr},size(activity));
    end
    
    %batch normalization
    if cnet.batchNorm(lyr)
        activity = Codat.NN.batchNormalization(activity, cnet.bn_beta{lyr}, ...
            cnet.bn_gamma{lyr}, cnet.bn_muInf{lyr}, cnet.bn_sig2Inf{lyr}, false);
    end
    
    %non_linearity
    activity = cnet.nonLinearity{lyr}(activity);
    
    %save activity shortcut
    idx = find(cnet.shortcut == lyr);
    if ~isempty(idx)
        activityShortcut{idx} = activity;
    end
end

prediction = activity;


end