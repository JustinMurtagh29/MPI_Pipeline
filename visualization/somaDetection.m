%% Read data mag 8
addpath(genpath('C:\code\'));
raw = readKnossosRoi('I:\CortexConnectomics\shared\cortex\2012-09-28_ex145_07x2\8\', '2012-09-28_ex145_07x2_mag8', [401 700; 307 606; 1 300]);

%% Remove blood vessels
bw = raw > 155 | raw < 70;
bw =  bwareaopen(bw, 1e4, 26);
bw = imclose(bw, ones(10,10,10));
bw = imfill(bw, 'holes');
raw(bw) = 122;

%% Prepare data
raw = single(raw);
raw = raw - mean(raw(:));
raw = raw ./ std(raw(:));

%% band pass filter (lowpass more important)
raw_filtered = gaussianbpf(raw,1,10); % detect large objects, filter out high frequencies
raw_filtered_2 = gaussianbpf(raw, 100, 1000); % maybe start using band pass filter in higher frequencies?, because they will be less in nuclei

%% Threshold
bw = raw_filtered > .3;

%% Exclude small objects
props = regionprops(bw, {'Area' 'PixelIdxList'});
excludeIdx = find([props(:).Area] < 1e5);
for i=1:length(excludeIdx)
    bw(props(excludeIdx(i)).PixelIdxList) = 0;
end

%% Make movie for visualization
makeClassMovie(bw, raw);
