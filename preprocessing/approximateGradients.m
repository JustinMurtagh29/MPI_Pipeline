function raw_smoothed = approximateGradients(vesselsMasked, bbox, filterSize)

raw = readKnossosRoi(vesselsMasked.root, vesselsMasked.prefix, bbox);

% Mean downsampling
display('Mean downsampling');
tic;
raw_mean = nlfilter3(raw, @mean, filterSize);
toc;

% Gradient calculation
display('Approximating gradients');
tic;
raw_temp = raw_mean;
% Mask out somata, apical dendrites and myelin (vessel already masked)
raw_temp(raw_temp < 100 | raw_temp > 130) = 121;
raw_smoothed = smooth3(raw_temp, 'gaussian', [5 5 5], 1.2);
toc;

end
