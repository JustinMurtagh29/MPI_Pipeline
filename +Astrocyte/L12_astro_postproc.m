astro_annot = load('~/GABA/astrocyte/predictions/unet_aug/Dropout02_val_v4_prediction.mat');
astro_vol_orig = int32(astro_annot.pred);

[x,y,z] = ndgrid(-4:4);
se = strel(sqrt(x.^2 + y.^2 + z.^2) <=4);

astro_vol = imclose(astro_vol_orig,se); 
astro_vol = bwareaopen(astro_vol, 5000,26);

pred = astro_vol; prediction = astro_vol_orig;
save('~/GABA/astrocyte/predictions/unet_aug/val_v4_03_postproc.mat', 'pred', 'prediction')
