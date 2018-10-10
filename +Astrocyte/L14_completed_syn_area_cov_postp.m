%{ 
author: Yagmur Yener
email: yagmur.yener.yy@gmail.com

run locally

astrocyte coverage synaptic cleft periphery vs synaptic cleft area
%}

% Load all synapses, and segments in small region

syn = load('~/GABA/astrocyte/synapses/syn.mat');
syn = syn.syn;

seg = load('~/GABA/astrocyte/synapses/segLarge1.mat');
seg = seg.seg;

astro_annot = load('~/GABA/astrocyte/predictions/unet_aug/v4_large1_postproc.mat');
astro_vol = int32(astro_annot.pred);

asiT = load('~/GABA/astrocyte/synapses/asiT.mat');
asiT = asiT.asiT;

preSynAxon = load('~/GABA/astrocyte/synapses/synID_axonSeg.mat');
preSynAxon = preSynAxon.syn_axon_segs;

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

%merge segments of pre and post
synSegments = cellfun(@vertcat, syn.synapses.presynId, syn.synapses.postsynId, 'UniformOutput', false);
%tell if the synapse id in the box given its total segments
lut_syn_idx = cellfun(@(synSegIds) all(lut_seg(synSegIds) == true), synSegments); %look-up table for syns

lut_syn_id = false([maxSynId,1]);
lut_syn_id(syn_areas(lut_syn_idx,1)) = true;

box_volume = numel(astro_vol) * 11.4/1000*11.4/1000*28/1000 ;
fprintf('Total of %d synapses in this %d um3 box.\n', sum(lut_syn_idx), round(box_volume))
lut_syn_idx = lut_syn_idx&(syn.synapses.type=='PrimarySpine');
fprintf('Total of %d primary spine synapses in this %d um3 box.\n', sum(lut_syn_idx), round(box_volume))


%% Locate synapses in volume

syn_idx = find(lut_syn_idx);
lut_syn_int = zeros(size(lut_syn_id));

for l = 1:sum(lut_syn_idx) 
    
    s = syn_idx(l); % syn idx
    
    ind_pre = []; ind_post = [];
    for i = 1:length(preSynAxon{s})
        ind_pre = [ind_pre; find(seg == preSynAxon{s}(i))]; %syn.synapses.presynId
    end
    
    for i = 1:length(postSynCompleted{s})
        ind_post =[ind_post; find(seg == postSynCompleted{s}(i))]; %syn.synapses.postsynId
    end
   
    mask = zeros(size(seg));
    mask(ind_pre) = 1;
    mask(ind_post) = 2;
    mask_pre = mask==1;
    mask_post = mask==2;

    mask_overlap=(imdilate(mask_pre,se)) & (imdilate(mask_post,se));
    
    mask_overlap_pre_post_astro = (mask_pre | mask_post) + 2*mask_overlap;
    mask_overlap_pre_post_astro(mask_overlap_pre_post_astro>1) = 2;
%     mask_overlap_pre_post_astro = mask_overlap_pre_post_astro + 5*double(imdilate(astro_vol,se));
    mask_overlap_pre_post_astro = mask_overlap_pre_post_astro + 5*double(astro_vol);
    mask_overlap_pre_post_astro(mask_overlap_pre_post_astro>2) = 5;
    
    % Find periphery of the cleft
    shiftedM = circshift(mask_overlap_pre_post_astro, [1,1]); %
    synPeriphery = (mask_overlap_pre_post_astro-shiftedM)==-3 | (mask_overlap_pre_post_astro-shiftedM)==2;
    astroSynInterface = mask_overlap_pre_post_astro-shiftedM==-3;
    
    shiftedM2 = circshift(mask_overlap_pre_post_astro, [-1,-1]); %
    synPeriphery = synPeriphery | ((mask_overlap_pre_post_astro-shiftedM2)==-3 | (mask_overlap_pre_post_astro-shiftedM2)==2);
    astroSynInterface = astroSynInterface | ((mask_overlap_pre_post_astro-shiftedM2)==-3);
           

    lut_syn_int(s) =  sum(astroSynInterface(:)==1)/sum(synPeriphery(:)==1)*100;
    
   
end

toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% plot area vs coverage of synapses
x=syn_areas(logical(lut_syn_idx~=0),2);
y=lut_syn_int(logical(lut_syn_idx~=0));

c=-linspace(1,10,length(x));
figure;colormap(lines);
scatter(x, y, 36, c, 'filled','MarkerFaceAlpha',.5,'MarkerEdgeAlpha',1)
% plot(syn_areas(logical(lut_syn~=0),2), lut_syn_int(logical(lut_syn~=0)), '*')
xlabel('Synaptic Cleft Area (um2)'); ylabel('Astrocyte Coverage (%)')

idx = isnan(y);
coefficients = polyfit(x(~idx), y(~idx), 1);
xFit = linspace(min(x), max(x), 1000);
yFit = polyval(coefficients , xFit);
hold on;
plot(xFit, yFit, 'k-');
% ax = gca;
% ax.YGrid = 'on';
% ax.XMinorGrid = 'on';
grid on