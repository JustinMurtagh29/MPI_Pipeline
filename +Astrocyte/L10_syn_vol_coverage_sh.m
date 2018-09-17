%{ 
author: Yagmur Yener
email: yagmur.yener.yy@gmail.com

run locally

Astrocyte Coverage vs synaptic volume with completed synapses (with shAgglos)
%}

% Load all synapses, and segments in small region

syn = load('~/GABA/astrocyte/synapses/syn.mat');
syn = syn.syn;
postSynCompleted = load('~/GABA/astrocyte/synapses/postSynCompleted.mat');
postSynCompleted=postSynCompleted.postSynCompleted;

seg = load('~/GABA/astrocyte/synapses/segLarge1.mat');
seg = seg.seg;

astro_annot = load('~/GABA/astrocyte/predictions/unet_aug/v4_large1.mat');
astro_vol = astro_annot.pred;

% Create a look-up table with all possible segment ids with segments in the
% loaded region are True

maxSegId = 15030572; %maximum possible segment ID

lut_seg = false(maxSegId, 1); %initialize as false
lut_seg(setdiff(seg, 0)) = true; %sets sed ids true if unique seg nonzero

% Find synapse IDs that are in this seg box (consists of segments)

%merge segments of pre and post
synSegments = cellfun(@vertcat, syn.synapses.presynId, postSynCompleted, 'UniformOutput', false);
%tell if the synapse id in the box given its total segments
lut_syn = cellfun(@(synSegIds) all(lut_seg(synSegIds) == true), synSegments); %look-up table for syns

box_volume = numel(astro_vol) * 11.4/1000*11.4/1000*28/1000 ;
fprintf('Total of %d synapses in this %d um3 box.\n', sum(lut_syn), round(box_volume))
lut_syn = lut_syn&syn.isSpineSyn;
fprintf('Total of %d primary spine synapses in this %d um3 box.\n', sum(lut_syn), round(box_volume))
%% Locate synapses in volume

synVolume = zeros(size(seg));
synVolume_l = synVolume;
mask = zeros(size(seg));

syn_idx = find(lut_syn);

for l = 1:sum(lut_syn) 
    
    s = syn_idx(l); % syn idx
    
    ind_pre = []; ind_post = [];
    for i = 1:length(syn.synapses.presynId{s})
        ind_pre = [ind_pre; find(seg == syn.synapses.presynId{s}(i))];
    end
    
    for i = 1:length(postSynCompleted{s})
        ind_post =[ind_post; find(seg == postSynCompleted{s}(i))];
    end
    
    % mark the location as primary spine
    mask(ind_pre) = 1;
    mask(ind_post) = 2;
    
    % mark the location with synapse idx
    synVolume(ind_pre) = s;
    synVolume(ind_post) = s;
    
    % just for simplifying the synapse ids a bit.
    synVolume_l(ind_pre) = l;
    synVolume_l(ind_post) = l;
    
end


%% Dilation 
% merging the segments of a synapse together.
% in synVolume, pre and post were also merged
% in 

% 3D object as dilation brush
[x,y,z] = ndgrid(-3:3);
se = strel(sqrt(x.^2 + y.^2 + z.^2) <=3);

% dilating-eroding the segments per synapse
synVolume_d = imdilate(synVolume_l,se); %note to future: use synVolume here
synVolume_de = imerode(synVolume_d, se);

% dilating-eroding segments of pre and post
mask_d = imdilate(mask,se);
mask_de = imerode(mask_d, se);


%  Separate pre and post synapses

mask_pre = mask_de; mask_post=mask_de/2;
mask_pre(mask_de==2) = 0;
mask_post(mask_de==1) = 0;
 

%% number of neighboring astrocyte pixels per synapse

% combine syn and astro in one mask (1:synapse, 5:astro)
mask_syn_astro = double(logical(mask_d)) + 5*double(astro_vol);
mask_syn_astro(mask_syn_astro>2) = 5;

% detect transition from 1 to 5 in xy plane
%make sure the border is within the area of the synapse segment so that one
%can get their overlap
shiftedM = circshift(mask_syn_astro, [1,1]); %
synPeriphery = (mask_syn_astro-shiftedM)==-4 | (mask_syn_astro-shiftedM)==1;
astroSynInterface = mask_syn_astro-shiftedM==-4;


shiftedM2 = circshift(mask_syn_astro, [-1,-1]); %
synPeriphery = synPeriphery | ((mask_syn_astro-shiftedM2)==-4 | (mask_syn_astro-shiftedM2)==1);
astroSynInterface = astroSynInterface | ((mask_syn_astro-shiftedM2)==-4);


%% count percent coverage of astrocytes of synapses


lut_syn_int = zeros(size(lut_syn));
synIds_d = unique(setdiff(synVolume_d, 0));

for i = 1:length(synIds_d)
    
    % a mask for one synapse id only
    test = zeros(size(synVolume));
    test(synVolume_d==synIds_d(i)) = 5;
    
    overlapInterfaceSyn = astroSynInterface+test;
    overlapPeripherySyn = synPeriphery+test;
    lut_syn_int(syn_idx(i)) =  sum(overlapInterfaceSyn(:)==6)/sum(overlapPeripherySyn(:)==6)*100;
    
end
%setdiff(lut_syn_int, 0)


%% find volume of synapses in um3

lut_syn_vol = zeros(size(lut_syn));
synIds_de = unique(setdiff(synVolume_de, 0));

for i = 1:length(synIds_de)

    lut_syn_vol(syn_idx(i)) = sum(synVolume_de(:)==synIds_de(i))  *11.4/1000*11.4/1000*28/1000 ;

end

%setdiff(lut_syn_vol, 0)

%% plot Volume vs coverage of synapses

plot(lut_syn_vol(logical(lut_syn~=0)), lut_syn_int(logical(lut_syn~=0)), '*')
xlabel('Synapse Volume (um3)'); ylabel('Astrocyte Coverage (%)')

