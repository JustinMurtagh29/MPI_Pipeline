function plotRawCube(root, prefix, bbox, outputFile)
    % Do not show figure (because no GUI output on cluster)
    figure('Visible', 'off');
    hold on;    
    bbox = repmat({bbox'},[6 1]);
    for i=1:6
        if mod(i,2)
            bbox{i}(i) = bbox{i}(i+1);
        else
            bbox{i}(i) = bbox{i}(i-1);
        end
        bbox{i} = bbox{i}';
        face{i} = readKnossosRoi(root,prefix,bbox{i});
        if i > 4
            face{i} = face{i}';
        end
        [X Y Z] = meshgrid(bbox{i}(1,1):bbox{i}(1,2),bbox{i}(2,1):bbox{i}(2,2),bbox{i}(3,1):bbox{i}(3,2));
        surf(squeeze(X),squeeze(Y),squeeze(Z),squeeze(face{i}),'facecolor','texture','edgecolor','none');
    end
    % Make it nicer looking xD
    colormap('gray');
    daspect([25 25 12]);
    view(3);
    camlight('headlight');
    lighting phong;
    axis off; axis tight;
    % One can control much more, e.g. Light color, strength, angle etc. -> google: Matlab Lightning Overview
    hgsave(gcf, outputFile, '-v7.3');
end

