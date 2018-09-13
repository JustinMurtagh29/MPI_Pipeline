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
lut_syn = cellfun(@(synSegIds) all(lut_seg(synSegIds) == true), synSegments); %look-up table for syns

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
    imagesc(synVolume_de(:,:,z)+double(astro_vol(:,:,z))*37, [0, 37]); colorbar off
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

mask_pre = mask_de; mask_post=mask_de/2;
mask_pre(mask_de==2) = 0;
mask_post(mask_de==1) = 0;

figure; colormap jet
for z = 1:72
    imagesc(synVolume_d(:,:,z)+double(astro_vol(:,:,z))*37, [0, 37]); colorbar
    pause(0.5)
end

%% Plot mask and astrocytes together

figure; colormap jet
for z = 1:72
    imagesc(mask_d(:,:,z)*12+double(astro_vol(:,:,z))*37, [0, 37]); colorbar
    pause(0.5)
end

%% number of neighboring astrocyte pixels per synapse


% combine syn and astro in one mask (1:synapse, 2:astro)
mask_syn_astro = double(logical(mask_d)) + 5*double(astro_vol);
mask_syn_astro(mask_syn_astro>2) = 5;

% detect transition from 1 to 2 in xy plane
shiftedM = circshift(mask_syn_astro, [1,1]);

% synPeriphery = abs(mask_syn_astro-shiftedM)==4 | abs(mask_syn_astro-shiftedM)==1;
% astroSynInterface = abs(mask_syn_astro-shiftedM)==4;

synPeriphery = (mask_syn_astro-shiftedM)==-4 | (mask_syn_astro-shiftedM)==-1;
astroSynInterface = mask_syn_astro-shiftedM==-4;


shiftedM2 = circshift(mask_syn_astro, [-1,-1]);
synPeriphery = synPeriphery | ((mask_syn_astro-shiftedM2)==-4 | (mask_syn_astro-shiftedM2)==1);
astroSynInterface = astroSynInterface | ((mask_syn_astro-shiftedM2)==-4);

%%
figure; colormap jet
for z = 1:72
    imagesc(abs(mask_syn_astro(:,:,z)-shiftedM(:,:,z))); colorbar
    pause(0.5)
end

%% count the number of pixels of astrocyte interface per synapse
% test1 = zeros(size(synVolume));
% test1(synVolume_d==14) = 5;
% figure; colormap jet
% imagesc(astroSynInterface(:,:,72)+ test1(:,:,72)); colorbar


lut_syn_int = zeros(size(lut_syn));
synIds_l = unique(setdiff(synVolume_d, 0));
synIds = unique(setdiff(synVolume, 0));

for i = 1:length(synIds_l)
    test = zeros(size(synVolume));
    test(synVolume_d==synIds_l(i)) = 5;
    overlapInterfaceSyn = astroSynInterface+test;
    overlapPeripherySyn = synPeriphery+test;
    lut_syn_int(synIds(i)) =  sum(overlapInterfaceSyn(:)==6)/sum(overlapPeripherySyn(:)==6);
    
end
setdiff(lut_syn_int, 0)
%%
 
figure; colormap jet
for z = 1:72
    imagesc(overlapPeripherySyn(:,:,z)); colorbar
    pause(0.5)
end