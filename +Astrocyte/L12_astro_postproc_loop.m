filenames1 = { "~/GABA/astrocyte/predictions/unet_aug/LR0001_Dropout02_Weight10_val2_prediction.mat", ...
               "~/GABA/astrocyte/predictions/unet_aug/LR0001_Dropout02_Weight10_val_prediction.mat", ...
               "~/GABA/astrocyte/predictions/unet_aug/LR0001_Dropout02_Weight10_val3_prediction.mat"};
filenames2 = { "~/GABA/astrocyte/predictions/unet_aug/LR0001_Dropout02_Weight10_val4_prediction.mat", ...
              "~/GABA/astrocyte/predictions/unet_aug/LR0001_Dropout02_Weight10_val5_prediction.mat"};
       
filenamesPost ={'~/GABA/astrocyte/predictions/unet_aug/val2_03_postproc.mat', ...
                '~/GABA/astrocyte/predictions/unet_aug/val_03_postproc.mat', ...
                '~/GABA/astrocyte/predictions/unet_aug/val3_03_postproc.mat', ...
                '~/GABA/astrocyte/predictions/unet_aug/val4_03_postproc.mat', ...
                '~/GABA/astrocyte/predictions/unet_aug/val5_03_postproc.mat'};
          
[x,y,z] = ndgrid(-4:4);
se = strel(sqrt(x.^2 + y.^2 + z.^2) <=4);
for i = 1:5
    
    if i <= 3
        predictions = load(filenames1{i});          
        astro_vol_orig = double(squeeze(predictions.pred) < 0);
    else
        predictions = load(filenames2{i-3});
        astro_vol_orig = double(squeeze(predictions.pred));
    end
    
    astro_vol = imclose(astro_vol_orig,se);
    astro_vol = bwareaopen(astro_vol, 5000,26);
    
    pred = astro_vol; prediction = astro_vol_orig;
    save(filenamesPost{i}, 'pred', 'prediction')
    
end




