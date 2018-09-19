astro_annot = load('~/GABA/astrocyte/predictions/unet_aug/v4_sw_val.mat');
astro_vol_orig = int32(astro_annot.pred);

[x,y,z] = ndgrid(-4:4);
se = strel(sqrt(x.^2 + y.^2 + z.^2) <=4);

astro_vol = imclose(astro_vol_orig,se); 
astro_vol = bwareaopen(astro_vol, 5000,26);

pred = astro_vol; prediction = astro_vol_orig;
save('~/GABA/astrocyte/predictions/unet_aug/v4_sw_val_postproc.mat', 'pred', 'prediction')
