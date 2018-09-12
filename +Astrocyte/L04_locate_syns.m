%{ 
author: Yagmur Yener
email: yagmur.yener.yy@gmail.com

run locally

Locating pre- and post- synapses in the astrocyte annotated segmented
volume
%}

% Load all synapses, and segments in small region

syn = load('~/GABA/astrocyte/synapses/syn.mat');
syn = syn.syn;

seg = load('~/GABA/astrocyte/synapses/seg.mat');
seg = seg.seg;

astro_annot = load('~/GABA/astrocyte/predictions/unet_aug/v4_val.mat');
astro_vol = astro_annot.pred;

% Create a look-up table with all possible segment ids with segments in the
% loaded region are True

maxSegId = 15030572; %maximum possible segment ID

lut_seg = false(maxSegId, 1); %initialize as false
lut_seg(setdiff(seg, 0)) = true; %sets sed ids true if unique seg nonzero

% Find synapse IDs that are in this seg box (consists of segments)

%merge segments of pre and post
synSegments = cellfun(@vertcat, syn.synapses.presynId, syn.synapses.postsynId, 'UniformOutput', false);
%tell if the synapse id in the box given its total segments
lut_syn = cellfun(@(synSegIds) any(lut_seg(synSegIds) == true), synSegments); %look-up table for syns

fprintf('Total of %d synapses in this box.\n', sum(lut_syn))
lut_syn = lut_syn&syn.isSpineSyn;
fprintf('Total of %d primary spine synapses in this box.\n', sum(lut_syn))
%% Locate synapses in volume

synVolume = zeros(size(seg));
synVolume_l = synVolume;
mask = zeros(size(seg));

syn_idx = find(lut_syn);

for l = 1:sum(lut_syn) %1:28
    
    s = syn_idx(l); % syn idx
    
    ind_pre = []; ind_post = [];
    for i = 1:length(syn.synapses.presynId{s})
        ind_pre = [ind_pre; find(seg == syn.synapses.presynId{s}(i))];
    end
    
    for i = 1:length(syn.synapses.postsynId{s})
        ind_post =[ind_post; find(seg == syn.synapses.postsynId{s}(i))];
    end
    
    % mark the location as primary spine
    mask(ind_pre) = 1;
    mask(ind_post) = 2;
    
    % mark the location with synapse idx
    synVolume(ind_pre) = s;
    synVolume(ind_post) = s;
    
    % just for simplifying the synapse ids a bit.
    synVolume_l(ind_pre) = l+5;
    synVolume_l(ind_post) = l+5;
    
end


%% Plot synapses and astrocytes together
% dark red is astrocytes

figure; colormap jet
for z = 1:72
    imagesc(synVolume_l(:,:,z)+double(astro_vol(:,:,z))*37, [0, 37]); colorbar
    pause(0.5)
end

%% Dilation 
% merging the segments of a synapse together.
% in synVolume, pre and post were also merged
% in 

% 3D object for dilation brush
[x,y,z] = ndgrid(-3:3);
se = strel(sqrt(x.^2 + y.^2 + z.^2) <=3);

% dilating-eroding the segments per synapse
synVolume_d = imdilate(synVolume_l,se);
synVolume_de = imerode(synVolume_d, se);
figure;
imshowpair(synVolume_l(:,:,20),synVolume_de(:,:,20),'montage')

% dilating-eroding segments of pre and post
mask_d = imdilate(mask,se);
mask_de = imerode(mask_d, se);
figure;
imshowpair(mask(:,:,20),mask_de(:,:,20),'montage')
%%  Plot either pre or post dilated synapses and astrocytes together

mask_pre = mask_d; mask_post=mask_d/2;
mask_pre(mask_d==2) = 0;
mask_post(mask_d==1) = 0;

figure; colormap jet
for z = 1:72
    imagesc(synVolume_d(:,:,z).*mask_post(:,:,z)+double(astro_vol(:,:,z))*37, [0, 37]); colorbar
    pause(0.5)
end

%% Plot mask and astrocytes together

figure; colormap jet
for z = 1:72
    imagesc(mask_d(:,:,z)*12+double(astro_vol(:,:,z))*37, [0, 37]); colorbar
    pause(0.5)
end

%% Define a measure for astrocyte coverage of synapses


