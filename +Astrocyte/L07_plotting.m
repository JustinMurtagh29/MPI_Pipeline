%{ 
author: Yagmur Yener
email: yagmur.yener.yy@gmail.com

run locally

plotting script
%}

%%
figure; colormap jet
z=5;
imagesc(double(astro_vol(:,:,z))*37, [0, 37]); colorbar
%% Dilated synapse volumes with astrocytes plotting 
figure; colormap jet
for z = 1:72
    imagesc(synVolume_d(:,:,z)*3+double(astro_vol(:,:,z))*37, [0, 37]); colorbar
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
% dark red: astrocyte
% light blue: synapse
% pink: boundary


figure; colormap jet
for z = 1:72
    imagesc(abs(mask_syn_astro(:,:,z)-shiftedM(:,:,z))); colorbar;
    pause(0.5)
end

%% Plot boundaries, astrocyte, synapse, whole synapse, overlap
vol = double(logical(synVolume_d))*5;
%light red: astro interface
%dark+light red: all interface

figure; colormap jet
for z = 1:72
    imagesc(vol(:,:,z) - astroSynInterface(:,:,z) + 2*synPeriphery(:,:,z), [0,7]); colorbar;
    pause(0.5)
end

%% plot boundaries and overlap with segments at the same time
vol = double(logical(synVolume_d))*5;

figure(1); colormap jet
figure(2); colormap jet
for z = 1:72
    figure(1);
    imagesc(vol(:,:,z) - astroSynInterface(:,:,z) + 2*synPeriphery(:,:,z), [0,7]); colorbar;
    figure(2);
    imagesc(abs(mask_syn_astro(:,:,z)-shiftedM(:,:,z))); colorbar;
    pause(0.5)
end
%% an example of overlappring syn and astro
z=167;
figure; colormap jet
imagesc(synVolume_de(:,:,z)+double(astro_vol(:,:,z))*37, [0, 37]); colorbar

%% plot Volume vs coverage of synapses
figure
plot(lut_syn_vol(logical(lut_syn~=0)), lut_syn_int(logical(lut_syn~=0)), '*')
xlabel('Synapse Volume (um3)'); ylabel('Astrocyte Coverage (%)')

%% look at the outliers
%light red: astro interface
%dark+light red: all this syn interface
%dark&light blue: syn&astro interface of other synapses (dark=asto)

vals = setdiff(lut_syn_int(:),0);
outliar = vals(1)
id = find(syn_idx==find(lut_syn_int == outliar))+5;

% a mask for one synapse id only
test = zeros(size(synVolume));
test(synVolume_d==id) = 5;

overlapInterfaceSyn = astroSynInterface+test;
overlapPeripherySyn = synPeriphery+test;

[~,~,Z]=ind2sub(size(synVolume), find(synVolume_d==id));
z_sorted = unique(Z);
figure; colormap jet
for i = 1:numel(z_sorted)
    z=z_sorted(i);
    imagesc(test(:,:,z) -astroSynInterface(:,:,z) + 2*synPeriphery(:,:,z), [0, 7]); colorbar;
    pause(0.5)
end

%% Plot boundaries, astrocyte, synapse, whole synapse, overlap
vol = synVolume_d/max(synVolume_d(:))*5;
%light red: astro interface
%dark+light red: all interface

figure; colormap jet
for z = 1:72
    imagesc(vol(:,:,z) + (- astroSynInterface(:,:,z) + 2*synPeriphery(:,:,z)), [0,7]); colorbar;
    pause(0.5)
end


%% 3D scatter plot surface vs volume vs coverage

figure
scatter3(lut_syn_vol(logical(lut_syn~=0)), syn_areas(logical(lut_syn~=0),2), lut_syn_int(logical(lut_syn~=0)), 'filled')
xlabel('Synapse Volume (um3)'); ylabel('Synaptic Cleft Area (um2)'); zlabel('Astrocyte Coverage (%)')

%% the removed marginal segments
figure;imshow(seg_mask(:,:,100))

%% copy two figures into one subfigure 

figure; colormap(lines)
h(1)=subplot(1,2,1);
title('Without Postprocessing')
xlabel('Synaptic Cleft Area (um2)'); ylabel('Astrocyte Coverage (%)')
h(2)=subplot(1,2,2);
title('With Postprocessing')
xlabel('Synaptic Cleft Area (um2)'); ylabel('Astrocyte Coverage (%)')
copyobj(allchild(get(figure(1),'Children')),h(1));
copyobj(allchild(get(figure(2),'CurrentAxes')),h(2));
grid(h(1), 'on')
grid(h(2), 'on')
