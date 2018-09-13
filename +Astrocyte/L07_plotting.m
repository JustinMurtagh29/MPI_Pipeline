%{ 
author: Yagmur Yener
email: yagmur.yener.yy@gmail.com

run locally

plotting script
%}


%% Dilated eroded synapse volumes plotting 
figure; colormap jet
for z = 1:72
    imagesc(synVolume_l(:,:,z)*3+double(astro_vol(:,:,z))*50, [0, 50]); colorbar off
    pause(0.5)
end

%% Before after dilation

figure; %Synapses
imshowpair(synVolume_l(:,:,20),synVolume_de(:,:,20),'montage')
figure; %Pre Post mask
imshowpair(mask(:,:,20),mask_de(:,:,20),'montage')

%% Plot synapses and astrocytes together
figure; colormap jet
for z = 1:72
    imagesc(synVolume_d(:,:,z)+double(astro_vol(:,:,z))*37, [0, 37]); colorbar
    pause(0.5)
end

%% Plot pre-post-mask and astrocytes together

figure; colormap jet
for z = 1:72
    imagesc(mask_d(:,:,z)*12+double(astro_vol(:,:,z))*37, [0, 37]); colorbar
    pause(0.5)
end

%% Plot the boundaries astrocyte, synapse, overlap

figure; colormap jet
for z = 1:72
    imagesc(abs(mask_syn_astro(:,:,z)-shiftedM(:,:,z))); colorbar;
    pause(0.5)
end

%% Plot boundaries, astrocyte, synapse, whole synapse, overlap
vol = double(logical(synVolume_d))*5;
% 2 astroSynInterface 
% 1 synPeriphery

figure; colormap jet
for z = 1:72
    imagesc(vol(:,:,z) + astroSynInterface(:,:,z) + synPeriphery(:,:,z), [0,7]); colorbar;
    pause(0.5)
end

%% plot boundaries and overlap with segments at the same time

figure(1); figure(2); colormap jet
for z = 1:72
    figure(1);
    imagesc(vol(:,:,z) + astroSynInterface(:,:,z) + synPeriphery(:,:,z), [0,7]); colorbar;
    figure(2);
    imagesc(abs(mask_syn_astro(:,:,z)-shiftedM(:,:,z))); colorbar;
    pause(0.5)
end

%% plot Volume vs coverage of synapses

plot(lut_syn_vol(logical(lut_syn~=0)), lut_syn_int(logical(lut_syn~=0)), '*')
xlabel('Synapse Volume (um3)'); ylabel('Astrocyte Coverage (%)')


