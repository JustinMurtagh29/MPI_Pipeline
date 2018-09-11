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
%%

synVolume = zeros(size(seg));
synVolume_l = synVolume;

mask_pre = zeros(size(seg)); mask_post = zeros(size(seg)); 
mask = zeros(size(seg));

syn_idx = find(lut_syn);

se = offsetstrel('ball',3,3);


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
    mask_pre(ind_pre) = 1;
    mask_post(ind_post) = 1;
    % dilate the masks
    mask_pre = imdilate(mask_pre,se);
    mask_post = imdilate(mask_post,se);
       
    % mark the location with synapse idx
    synVolume(logical(mask_pre)) = s;
    synVolume(logical(mask_post)) = s;
    
    % just for simplifying the synapse ids a bit.
    synVolume_l(logical(mask_pre)) = l;
    synVolume_l(logical(mask_post)) = l;
    
    % have a mask identifying where the pre and post synapses are (1:pre, 2:post, 0:bg)
    mask = mask + mask_pre + 2*mask_post; %assuming non-overlapping segments
    
    %reset the masks
    mask_pre = zeros(size(seg)); mask_post = zeros(size(seg)); 
    
end

%% Plot synapses and astrocytes together
% dark red is astrocytes

figure; colormap jet
for z = 1:72
    imagesc(synVolume_l(:,:,z)+5+double(astro_vol(:,:,z))*37, [0, 37]); colorbar
    pause(0.5)
end


