function t = thresholdsForFocusedAnnotation( rp, thresh, reqP, reqR )
%THRESHOLDSFORFOCUSEDANNOTATION Determine threshold for focused annotation
% when a given precision and recall rate is required.
% INPUT rp: [Nx2] double
%           Recall-recision values on the test set.
%       thresh: [Nx1] double
%           Thresholds respective rp values.
%       reqP: double
%           Required precision.
%       reqR: double
%           Required recall.
% OUTPUT t: [1x2] double
%           Lower and the upper threshold for focused annotation.
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

lowerThresIdx = find(rp(:,1) >= reqR, 1, 'last');
upperThresIdx = find(rp(:,2) >= reqP, 1, 'first');
t(1) = thresh(lowerThresIdx);
t(2) = thresh(upperThresIdx);

end

