load('D:\sync\activeTraining\gradientData.mat')

%% Visualize blood vessel masked and mean downsampled data
figure('Position', [1 1 1600 1124]);
colors = colormap(hot(128));
for i=1:size(raw_mean,3)
   subplot(4,4,[2 3 4 6 7 8 10 11 12]);
   imagesc(raw_mean(:,:,i)');
   caxis([60 160]);
   axis equal; axis off;
   subplot(4,4,[1 5 9]);
   hold on;
   plot(mean(raw_mean(:,:,i),1), 'Color', colors(i,:));
   view([90 90]);
   axis off;
   xlim([0 112]); ylim([110 130]);
   subplot(4,4,[14 15 16]);
   hold on;
   plot(mean(raw_mean(:,:,i),2), 'Color', colors(i,:));
   axis off;
   xlim([0 160]); ylim([110 130]);
   colormap(colors);
   drawnow;
   pause(.5);
end
%% Plot gradients along all directions

along{1} = squeeze(mean(mean(raw_mean,2),3));
along{2} = squeeze(mean(mean(raw_mean,1),3));
along{3} = squeeze(mean(mean(raw_mean,1),2));
figure('Position', [1 41 1600 1084]);
for i=1:3
    subplot(3,1,i);
    plot(along{i});
end

%% Plot single layers
colors = colormap(hot(128));
for i=1:size(raw_mean,3)
    subplot(2,1,1);
    hold on;
    plot(mean(raw_mean(:,:,i),1), 'Color', colors(i,:));
    subplot(2,1,2);
    hold on;
    plot(mean(raw_mean(:,:,i),2), 'Color', colors(i,:));
end

%% Look at weird peaks vs draughts
figure('Position', [1 41 1600 1084]);
yPos = [21 26 29 36 42];
for i=1:length(yPos)
    subplot(1,length(yPos),i);
    imagesc(squeeze(raw_smoothed(:,yPos(i),:)));
    caxis([60 160]);
    axis equal; axis off;
end

%% Smooth gradients
raw_temp = raw_mean;
raw_temp(raw_temp < 100 | raw_temp > 130) = 121; % Mask out somata, apical dendrites and myelin
raw_smoothed = smooth3(raw_temp, 'gaussian', [5 5 5], 1.2);

%% Visualize blood vessel masked and mean downsampled data
close all;
figure('Position', [1 1 1600 1124]);
colors = colormap(hot(128));
for i=1:size(raw_temp,3)
   subplot(4,4,[2 3 4 6 7 8 10 11 12]);
   imagesc(raw_temp(:,:,i)');
   caxis([100 130]);
   axis equal; axis off;
   subplot(4,4,[1 5 9]);
   hold on;
   plot(mean(raw_temp(:,:,i),1), 'Color', colors(i,:));
   view([90 90]);
   axis off;
   xlim([0 112]); ylim([110 130]);
   subplot(4,4,[14 15 16]);
   hold on;
   plot(mean(raw_temp(:,:,i),2), 'Color', colors(i,:));
   axis off;
   xlim([0 160]); ylim([110 130]);
   colormap(colors);
   drawnow;
   pause(.5);
end

%% Visualize blood vessel masked and mean downsampled data
close all;
figure('Position', [1 1 1600 1124]);
colors = colormap(hot(128));
for i=1:size(raw_smoothed,3)
   subplot(4,4,[2 3 4 6 7 8 10 11 12]);
   imagesc(raw_smoothed(:,:,i)');
   caxis([100 130]);
   axis equal; axis off;
   subplot(4,4,[1 5 9]);
   hold on;
   plot(mean(raw_smoothed(:,:,i),1), 'Color', colors(i,:));
   view([90 90]);
   axis off;
   xlim([0 112]); ylim([110 130]);
   subplot(4,4,[14 15 16]);
   hold on;
   plot(mean(raw_smoothed(:,:,i),2), 'Color', colors(i,:));
   axis off;
   xlim([0 160]); ylim([110 130]);
   colormap(colors);
   drawnow;
   pause(.5);
end

%% Plot single layers
colors = colormap(hot(128));
for i=1:size(raw_smoothed,3)
    subplot(2,1,1);
    hold on;
    plot(mean(raw_smoothed(:,:,i),1), 'Color', colors(i,:));
    subplot(2,1,2);
    hold on;
    plot(mean(raw_smoothed(:,:,i),2), 'Color', colors(i,:));
end

%% Plot gradients along all directions

along{1} = squeeze(mean(mean(raw_smoothed,2),3));
along{2} = squeeze(mean(mean(raw_smoothed,1),3));
along{3} = squeeze(mean(mean(raw_smoothed,1),2));
figure('Position', [1 41 1600 1084]);
for i=1:3
    subplot(3,1,i);
    plot(along{i});
end

%% TESTSUITE for smoothing param
%% Plot single layers
figure('Position', [1 41 1600 1084]);
colors = colormap(hot(128));
for i=1:size(raw_mean,3)
    subplot(2,1,1);
    hold on;
    plot(mean(raw_mean(:,:,i),2), 'Color', colors(i,:));
    subplot(2,1,2);
    hold on;
    plot(mean(raw_mean(:,:,i),1), 'Color', colors(i,:));
end
raw_smoothed = smooth3(raw_mean, 'gaussian', [7 7 7], 3);
figure('Position', [1601 1 1920 1124]);
for i=1:size(raw_smoothed,3)
    subplot(2,1,1);
    hold on;
    plot(mean(raw_smoothed(:,:,i),2), 'Color', colors(i,:));
    subplot(2,1,2);
    hold on;
    plot(mean(raw_smoothed(:,:,i),1), 'Color', colors(i,:));
end

