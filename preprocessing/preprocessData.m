% Reading raw data from HDF5 file
display('Loading data from HDF5 file');
tic;
raw = h5read('/gaba/scratch/mberning/tempData07x2.h5', '/raw');
toc;

display('Thresholding & Filling:');
tic;
% generate logical arrrays
bloodVesselMask = zeros(size(raw), 'uint8');
toc;

tic;
fprintf('\nLayer Counter (max 3422): ');
for i=1:size(raw,3)
    % Nice counter, i Like
    if i>1
		for j=0:log10(i-1)
			fprintf('\b');
		end
	end
	fprintf('%d', i);

	temp = bwareaopen(raw(:,:,i) > 162 | raw(:,:,i) < 50, 1000, 4);
	temp = imclose(temp, ones(5,5));	
	% Do we need to change location of hole drilling? -> test
    temp = padarray(temp, [1 1], 1);
    idx = 1;
	while temp(1800,idx) == 1
		temp(1800,idx) = 0;
		idx = idx + 1;
	end
	temp = imfill(temp, 'holes');
	temp([1 end], :) = [];
	temp(:, [1 end]) = [];
	bloodVesselMask(:,:,i) = temp;
end
fprintf('\n');
clear temp idx;
toc;

display('Saving blood vessel mask to HDF 5 file:');
tic;
h5write('/gaba/scratch/mberning/tempData07x2.h5', '/mask', bloodVesselMask);
toc;

display('Masking the raw data:');
tic;
fprintf('\nLayer Counter (max 3422): ')
for i=1:size(raw,3);
	if i>1
		for j=0:log10(i-1)
			fprintf('\b');
		end
	end
	fprintf('%d', i);
	temp = raw(:,:,i);
	temp(b(:,:,i)>0) = 121;
	raw(:,:,i) = temp;
end
fprintf('\n')
toc;

display('Saving masked raw data to HDF 5 file:');
tic;
raw = uint8(raw);
h5write('/gaba/scratch/mberning/tempData07x2.h5', '/rawmasked', raw);
toc;

% Mean downsampling
filterSize = [64 64 29];
display('Mean downsampling');
tic;
raw_mean = nlfilter3(raw, @mean, filterSize);
toc;

% Gradient calculation
display('Approximating gradients');
tic;
raw_temp = raw_mean;
raw_temp(raw_temp < 100 | raw_temp > 130) = 121; % Mask out somata, apical dendrites and myelin
raw_smoothed = smooth3(raw_temp, 'gaussian', [5 5 5], 1.2);
toc

display('Gradient correction');
tic;
fprintf('\nLayer Counter (max 3422): ')
for i=1:size(raw,3);
	if i>1
		for j=0:log10(i-1)
			fprintf('\b');
		end
	end
	fprintf('%d', i);
	z = ceil(i/filterSize(3));
	planeCorrectionImage = 1./(interp(mean(raw_smoothed(:,:,z),2),filterSize(1))*interp(mean(raw_smoothed(:,:,z),1),filterSize(2))./121^2);
	temp = raw(:,:,i);
	temp = planeCorrectionImage.*temp;
	temp(b(:,:,i)>0) = 121;
	raw(:,:,i) = temp;
end
fprintf('\n')
toc;

display('Saving gradient masked data to HDF 5 file:');
tic;
raw = uint8(raw);
h5create('/gaba/scratch/mberning/tempData07x2.h5', '/rawcorrected', size(raw), 'Datatype', 'uint8');
h5write('/gaba/scratch/mberning/tempData07x2.h5', '/rawcorrected', raw);
toc;

display('Writing to KNOSSOS hierachy');
tic;
writeKnossosRoi('/zdata/manuel/data/cortex/2012-09-28_ex145_07x2/mag1masked/', '2012-09-28_ex145_07x2_mag1', bbox(:,1)', raw);
toc;

