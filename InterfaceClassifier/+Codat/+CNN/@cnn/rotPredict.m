function prediction = rotPredict( cnet, input )
%ROTPREDICT Prediction for rotation invariant feature cnet.
% INPUT input: Input cube
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

cnet.isTraining = false;
prediction = zeros(size(input) - cnet.border,'like',input);
for rotIter = 1:4
    prediction = prediction + ...
                 rot90(cnet.predict( rot90(input,rotIter) ), 4-rotIter)./4;
end

end
