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

%% Astrocyte Volume Percentage per Slice
figure; hold on 
idx = [40:160];
perc1 = zeros(1,size(astro_vol,3));
for i = 1:size(astro_vol,3)
    perc1(i) = sum(sum(astroSynInterface(:,:,i)==1) / sum(synPeriphery(:,:,i)==1));
end
plot(idx, perc1(idx)/mean(perc1))

% perc = zeros(1,size(astro_vol,3));
% for i = 1:size(astro_vol,3)
%     perc(i) = sum(sum(synPeriphery(:,:,i)==1));
% end
% plot([40:160], perc([40:160])/max(perc))

perc = zeros(1,size(astro_vol,3));
for i = 1:size(astro_vol,3)
    perc(i) = sum(sum(astro_vol(:,:,i)==1))/ (size(astro_vol,1)*size(astro_vol,2)) *100;
end
plot(idx, perc(idx)/mean(perc))

legend('AstroBorder Perc', 'Astro Density', 'Location', 'southeast')
ylabel('Boundary')
xlabel('Z slice')
grid on

figure;
colors=[1:size(astro_vol,3)];
scatter(perc(idx), perc1(idx), 36,idx,'filled')
xlabel('Astro Density'); ylabel('AstroBorderPerc')
h = colorbar;
ylabel(h, 'Z slices')
%% Synapse alongation in z versus astro coverage

figure; hold on
for i = 1:length(synIds_d)
        
    % a mask for one synapse id only
    test = zeros(size(synVolumeOrig));
    test(synVolume_d==synIds_d(i)) = 5;
    
    plot(lut_syn_int(syn_idx(i))*(squeeze(sum(sum(test, 1),2))>1), '.')

%     overlapInterfaceSyn = astroSynInterface+test;
%     overlapPeripherySyn = synPeriphery+test;
%     lut_syn_int(syn_idx(i)) =  sum(overlapInterfaceSyn(:)==6)/sum(overlapPeripherySyn(:)==6)*100;
    
end

%% Look at one synapse pre and post astro coverage

vals = setdiff(lut_syn_int(:),0);
outliar = vals(1)
id = find(syn_idx==find(lut_syn_int == outliar))+5;

% a mask for one synapse id only
test = zeros(size(synVolume));
test(synVolume_d==id) = 1;

overlapInterfaceSyn = astroSynInterface+test;
overlapPeripherySyn = synPeriphery+test;

[~,~,Z]=ind2sub(size(synVolume), find(synVolume_d==id));
z_sorted = unique(Z);
figure(1); colormap jet
for i = 1:numel(z_sorted)
    z=z_sorted(i);
%     imagesc(test(:,:,z) -astroSynInterface(:,:,z) + 2*synPeriphery(:,:,z), [0, 7]); colorbar;
    imagesc(test(:,:,z).*mask_d(:,:,z)*12+double(astro_vol(:,:,z))*37, [0, 37]); colorbar

    pause(0.5)
end

%% RAW OVERLAY

raw = load('~/GABA/astrocyte/synapses/raw_val.mat');
raw = raw.raw;
E = raw(:,:,1);
red = cat( 3, ones(size(E)), zeros(size(E)), zeros(size(E)) );
blue = cat( 3, zeros(size(E)), zeros(size(E)), ones(size(E)) );
green = cat( 3, zeros(size(E)), ones(size(E)), zeros(size(E)) );
magenta = cat( 3, ones(size(E)), zeros(size(E)), 3*ones(size(E)) );
%% Transparent mask overlay on Raw
% figure; colormap gray
for z = 1:72
    a = imshow(raw(:,:,z), [], 'InitialMag', 'fit');
    hold on
    r = imshow(red);
    r.AlphaData = logical(astro_vol(:,:,z))*0.3;
    b = imshow(blue);
    b.AlphaData = mask_de(:,:,z)*0.15;
    hold off
    pause(0.5)
end

%% Raw overlay with boundaries

figure; colormap gray
for z = 1:72
    imagesc(raw(:,:,z))%, 'InitialMag', 'fit');
    hold on
    r = imshow(red);
    r.AlphaData = (abs(mask_syn_astro(:,:,z)-shiftedM(:,:,z))==5);
    b = imshow(blue);
    b.AlphaData = (abs(mask_syn_astro(:,:,z)-shiftedM(:,:,z))==1);
    m = imshow(magenta);
    m.AlphaData = (abs(mask_syn_astro(:,:,z)-shiftedM(:,:,z))==4);
    hold off
    pause(0.5)
end

%% Raw overlay with pre and post overlap
% [x,y,z] = ndgrid(-2:2);
% se = strel(sqrt(x.^2 + y.^2 + z.^2) <=2);
% % mask_overlap=((imdilate((mask_de==1).*mask_syn,se)) & (imdilate((mask_de==2).*mask_syn,se)));
% mask_overlap=(imdilate(mask_pre,se)) & (imdilate(mask_post,se));

figure(1); colormap gray
for z = 10:52
    a = imshow(raw(:,:,z), [], 'InitialMag', 'fit');
    hold on
    r = imshow(red);
    r.AlphaData = logical(astro_vol(:,:,z))*0.3;
    b = imshow(blue);
    b.AlphaData = (mask(:,:,z))*0.15;
    g = imshow(green);
    g.AlphaData = mask_overlap(:,:,z)*0.5;
    hold off
    pause(0.5)
end

%% mask_overlap_pre_post_astro

figure(); colormap gray
for z = 10:45
    a = imshow(raw(:,:,z), [], 'InitialMag', 'fit');
    hold on
    r = imshow(red);
    r.AlphaData = (mask_overlap_pre_post_astro(:,:,z)==5)*0.3;
    b = imshow(blue);
    b.AlphaData = (mask_overlap_pre_post_astro(:,:,z)==1)*0.15;
    g = imshow(green);
    g.AlphaData = (mask_overlap_pre_post_astro(:,:,z)==2)*0.5;
    hold off
    pause(0.5)
end