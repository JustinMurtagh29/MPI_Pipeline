close all; clc;
addpath('/home/mberning/code/auxiliaryMethods/');
cortexDataset = '/media/big2/data/2012-06-08_2012-06-04/mag1/';
saveDir = '/home/mberning/Desktop/prototyp1/';
startPoint = [1600 1300 510];

%%
bbox = [startPoint - [100 100 0]; startPoint + [100 100 99]];
bbox = bbox';
figure('position', [1601 2 1920 1120]);
for i=1:10
    raw = readKnossosRoi(cortexDataset, '2012-06-28_Cortex_mag1', bbox);
    for f=1:size(raw,3)
        imagesc(raw(:,:,f));
        colormap('gray');
        axis equal;
        axis off;        
        title([num2str(i) ':' num2str(f)]);
        if i ==1 && f == 1 
           [x,y] = ginput(1);
        elseif f == 1
            switch exitDir
                case 'o'
                    x = x2;
                    y = 201;
                case 'u'
                    x = x2;
                    y = 1;
                case 'l'
                    x = 201;
                    y = y2;
                case 'r'
                    x = 1;
                    y = y2;
                otherwise
                    x = x2;
                    y = y2;
            end
            hold on;
            plot(x,y, 'or', 'MarkerSize', 10, 'LineWidth', 4);
            hold off;
            pause;
        end
        pause(.3);
        if strcmp(get(gcf, 'CurrentCharacter'), 's')
            break
        end
    end
    exitDir = input('Where did it exit (o/u/r/l/n)?', 's');
    exitFrame = f;
    switch exitDir
        case 'o'
            bbox(1,:) = bbox(1,:) - [200 200];
            bbox(3,:) = bbox(3,:) + [exitFrame exitFrame];
        case 'u'
            bbox(1,:) = bbox(1,:) + [200 200];
            bbox(3,:) = bbox(3,:) + [exitFrame exitFrame];
        case 'l'
            bbox(2,:) = bbox(2,:) - [200 200];
            bbox(3,:) = bbox(3,:) + [exitFrame exitFrame];
        case 'r'
            bbox(2,:) = bbox(2,:) + [200 200];
            bbox(3,:) = bbox(3,:) + [exitFrame exitFrame];
        otherwise
            bbox(3,:) = bbox(3,:) + [exitFrame exitFrame];
    end
    if ~exist([saveDir num2str(i, '%2.2i') '/'], 'dir')
        mkdir([saveDir num2str(i, '%2.2i') '/']);
    end
    for f=1:exitFrame
        imagesc(raw(:,:,f));
        colormap('gray');
        axis equal;
        axis off;
        title([num2str(i) ':' num2str(f)]);
        if f == 1
            hold on;
            plot(x,y, 'or', 'MarkerSize', 10, 'LineWidth', 4);
            pause(.5);
            hold off;
        elseif f == exitFrame
            [x2,y2] = ginput(1);
            hold on;
            plot(x2,y2, 'or', 'MarkerSize', 10, 'LineWidth', 4);
            hold off;
        end
        saveas(gcf, [saveDir num2str(i, '%2.2i') '/image' num2str(f, '%3.3i') '.tif']);
    end
end