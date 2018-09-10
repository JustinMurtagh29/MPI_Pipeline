%{ 
author: Yagmur Yener
email: yagmur.yener.yy@gmail.com

run locally

Loading and looking at the file that was saved using S01_load_and_save
%}

%% Load the predictions from file

%location in the L4 dataset
box_offset = [4179, 4994, 2264];
box_shape = [178, 178, 72];
bbox = [box_offset' , (box_offset+box_shape)'];

%Region annotated by CNN (voxel-wise)
%pred: binary thresholded predictions 1 means astro
%prediction: predicted probabilities continuous between -1.7 and 1.7
%(scaled tanh output). -1 means astro
astro_annot = load('~/GABA/astrocyte/predictions/unet_aug/v4_val.mat');
astro_vol = astro_annot.pred;

%% Load synapse point from file

syn_points = load('~/GABA/astrocyte/synapses/syn_points.mat');
syn_points = syn_points.syn_points;

syn = load('~/GABA/astrocyte/synapses/syn.mat');
syn = syn.syn;
conn = load('~/GABA/astrocyte/synapses/conn.mat');
conn = conn.conn;
param = load('~/GABA/astrocyte/synapses/param.mat');
param = param.param;


%% Astrocyte volume to points

[x, y, z, ~] = ind2sub(size(astro_vol), find(astro_vol)); %linear to 3D indices
astro_points = [x, y, z] + box_offset; 

%% Retrieve astrocyte volume from points

test_astro_vol = zeros(box_shape);
test_astro_vol( sub2ind(box_shape, astro_points(:,1)-box_offset(1), astro_points(:,2)-box_offset(2), astro_points(:,3)-box_offset(3)) ) = 1;
isequal(test_astro_vol, astro_vol)

%% Retrieve synapse volume from points

%ind2sub whole

%cut with offset and shape




