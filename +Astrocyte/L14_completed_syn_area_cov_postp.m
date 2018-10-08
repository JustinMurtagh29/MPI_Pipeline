%{ 
author: Yagmur Yener
email: yagmur.yener.yy@gmail.com

run locally

astrocyte coverage synaptic cleft periphery vs synaptic cleft area
%}

% Load all synapses, and segments in small region

syn = load('~/GABA/astrocyte/synapses/syn.mat');
syn = syn.syn;

seg = load('~/GABA/astrocyte/synapses/seg.mat');
seg = seg.seg;

astro_annot = load('~/GABA/astrocyte/predictions/unet_aug/v4_val.mat');
astro_vol = int32(astro_annot.pred);

asiT = load('~/GABA/astrocyte/synapses/asiT.mat');
asiT = asiT.asiT;

preSynAxon = load('~/GABA/astrocyte/synapses/preSynAxon.mat');
preSynAxon = preSynAxon.preSynAxon;

postSynCompleted = load('~/GABA/astrocyte/synapses/postSynCompleted.mat');
postSynCompleted = postSynCompleted.postSynCompleted;

maxSegId = 15030572; %maximum possible segment ID

%%
tic

%Segments that are not leaving the cube (remove marginals)
seg_mask = uint32(seg>0);
seg_mask = imclearborder(seg_mask,6); %connectivity in 3D
seg=seg.*seg_mask;

% Smooth out the segments 
[x,y,z] = ndgrid(-3:3); % 3D object as dilation brush
se = strel(sqrt(x.^2 + y.^2 + z.^2) <=3);
seg = imclose(seg,se); 

% Segments that are in this volume
lut_seg = false(maxSegId, 1); %initialize as false
lut_seg(setdiff(seg, 0)) = true; %sets sed ids true if unique seg nonzero

% Synapse id to synapse area map
syn_areas = [asiT.id asiT.area];
syn_areas = sortrows(syn_areas);

% Remove some synapses
maxSynId = height(syn.synapses);
syn.synapses = syn.synapses(syn_areas(:,1) , :);
preSynAxon = preSynAxon(syn_areas(:,1));
postSynCompleted = postSynCompleted(syn_areas(:,1));
%%
%merge segments of pre and post
% synSegments = cellfun(@vertcat, syn.synapses.presynId, syn.synapses.postsynId, 'UniformOutput', false);
% %tell if the synapse id in the box given its total segments
% lut_syn_idx = cellfun(@(synSegIds) all(lut_seg(synSegIds) == true), synSegments); %look-up table for syns

%tell if the synapse id in the box given it's spine head
lut_syn_idx = cellfun(@(synSegIds) all(lut_seg(synSegIds) == true), syn.synapses.presynId);

lut_syn_id = false([maxSynId,1]);
lut_syn_id(syn_areas(lut_syn_idx,1)) = true;

box_volume = numel(astro_vol) * 11.4/1000*11.4/1000*28/1000 ;
fprintf('Total of %d synapses in this %d um3 box.\n', sum(lut_syn_idx), round(box_volume))
lut_syn_idx = lut_syn_idx&(syn.synapses.type=='PrimarySpine');
fprintf('Total of %d primary spine synapses in this %d um3 box.\n', sum(lut_syn_idx), round(box_volume))


%% Locate synapses in volume

synVolumeOrig = cell(size(seg)); %syn volume with original ids
synVolume = synVolumeOrig; %syn volume with simplified ids for images
mask = zeros(size(seg));

syn_idx = find(lut_syn_idx);

for l = 1:sum(lut_syn_idx) 
    
    s = syn_idx(l); % syn idx
    
    ind_pre = []; ind_post = [];
    for i = 1:length(preSynAxon{s})
        ind_pre = [ind_pre; find(seg == preSynAxon{s}(i))]; %syn.synapses.presynId
    end
    
    for i = 1:length(postSynCompleted{s})
        ind_post =[ind_post; find(seg == postSynCompleted{s}(i))]; %syn.synapses.postsynId
    end
   
    % mark the location as primary spine
    mask(ind_pre) = 1;
    mask(ind_post) = 2;
           
    % mark the location with synapse idx
    synVolumeOrig(ind_pre) = {s};
    synVolumeOrig(ind_post) = {s};
    
    % just for simplifying the synapse ids a bit.
    synVolume(ind_pre) = {[synVolume{ind_pre},l+5]};
    synVolume(ind_post) = {[synVolume{ind_post},l+5]};
    
   
end


%% Dilation 
% merging the segments of a synapse together.
% in synVolume, pre and post were also merged
% in 

% 3D object as dilation brush
[x,y,z] = ndgrid(-3:3);
se = strel(sqrt(x.^2 + y.^2 + z.^2) <=3);

% % dilating-eroding the segments per synapse
% synVolume_d = imdilate(synVolume,se); 
% synVolume_de = imerode(synVolume_d, se);

% dilating-eroding segments of pre and post
mask_d = imdilate(mask,se);
mask_de = imerode(mask_d, se);


%  Separate pre and post synapses
mask_pre = mask_de; mask_post=mask_de/2;
mask_pre(mask_de==2) = 0;
mask_post(mask_de==1) = 0;
mask_overlap=(imdilate(mask_pre,se)) & (imdilate(mask_post,se));
mask_overlap_pre_post = (imdilate(mask_pre,se)) + (imdilate(mask_post,se));

%% count percent coverage of astrocytes of synapses


lut_syn_int = zeros(size(lut_syn_id));
synIds_d = unique(setdiff(synVolume_d, 0));

for i = 1:length(synIds_d)
    
%     mask_syn = cellfun(@(vox) vox==synIds_d(i),synVolume,'UniformOutput',false);
        
    % combine syn and astro in one mask (1:synapse, 5:astro)
    mask_syn_astro = double(logical(synVolume_d==synIds_d(i))) + 5*double(astro_vol);
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
    
    % a mask for one synapse id only
    test = zeros(size(synVolumeOrig));
    test(synVolume_d==synIds_d(i)) = 5;
    
    overlapInterfaceSyn = astroSynInterface+test;
    overlapPeripherySyn = synPeriphery+test;
    lut_syn_int(syn_idx(i)) =  sum(overlapInterfaceSyn(:)==6)/sum(overlapPeripherySyn(:)==6)*100;
    
end
%setdiff(lut_syn_int, 0)
toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% number of neighboring astrocyte pixels to synapse clusters

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
%% synapse volume

lut_syn_vol = zeros(size(lut_syn_id));
synIds_de = unique(setdiff(synVolume_de, 0));

for i = 1:length(synIds_de)

    lut_syn_vol(syn_idx(i)) = sum(synVolume_de(:)==synIds_de(i))  *11.4/1000*11.4/1000*28/1000 ;

end

%% plot area vs coverage of synapses
x=syn_areas(logical(lut_syn_idx~=0),2);
y=lut_syn_int(logical(lut_syn_idx~=0));

c=-linspace(1,10,length(x));
figure;colormap(lines);
scatter(x, y, 36, c, 'filled','MarkerFaceAlpha',.5,'MarkerEdgeAlpha',1)
% plot(syn_areas(logical(lut_syn~=0),2), lut_syn_int(logical(lut_syn~=0)), '*')
xlabel('Synaptic Cleft Area (um2)'); ylabel('Astrocyte Coverage (%)')


coefficients = polyfit(x, y, 1);
xFit = linspace(min(x), max(x), 1000);
yFit = polyval(coefficients , xFit);
hold on;
plot(xFit, yFit, 'k-');
ax = gca;
ax.YGrid = 'on';
ax.XMinorGrid = 'on';
