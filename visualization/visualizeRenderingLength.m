function visualizeRenderingLength( v, id )

for j=1:length(v)
    for i=1:length(v(j).missions)
        % Rendering Length
        startFF(i) = v(j).missions(i).start.firstFrame;
        startLF(i) = v(j).missions(i).start.lastFrame;
        for k=1:length(v(j).missions(i).possibleEnds)
            endFF{i,k} = v(j).missions(i).possibleEnds(k).firstFrame;
            endLF{i,k} = v(j).missions(i).possibleEnds(k).lastFrame;
        end
        % Number possible ends
        start = v(j).missions(i).start.id;
        possibleEnds{i} = [v(j).missions(i).possibleEnds(:).id];
        nrPossibleEnds(i) = length(possibleEnds{i});
        probabilities{i} = [v(j).missions(i).possibleEnds(:).probability];
    end

    % Figure out probabilty cutoff
    for thres=1:15
        t(thres) = (thres*0.02);
        for i=1:length(v(j).missions)
            nrObjects(i,thres) = sum(probabilities{i} > t(thres));
        end
    end

    figure('Units', 'centimeters', 'Position', [0 0 29.7 21], ...
        'PaperType', 'A4', 'PaperUnits', 'centimeters', 'PaperPosition', [0 0 29.7 21], 'PaperOrientation', 'landscape');
    subplot(4,1,1);
    hist(startFF, 20);
    xlabel('distance along problem direction [nm]');
    ylabel('occurences');
    title('First frame start segment');
    subplot(4,1,2);
    hist(startLF, 20);
    xlabel('distance along problem direction [nm]');
    ylabel('occurences');
    title('Last frame start segment');
    subplot(4,1,3);
    hist([endFF{:}], 20);
    xlabel('distance along problem direction [nm]');
    ylabel('occurences');
    title('First frame end segment');
    subplot(4,1,4);
    hist([endLF{:}], 20);
    xlabel('distance along problem direction [nm]');
    ylabel('occurences');
    title('Last frame end segment');
    saveas(gcf, ['C:\Users\mberning\Desktop\problemInspector\' id v(j).name '000RenderingLength.pdf']);
    close all;

    figure('Units', 'centimeters', 'Position', [0 0 29.7 21], ...
        'PaperType', 'A4', 'PaperUnits', 'centimeters', 'PaperPosition', [0 0 29.7 21], 'PaperOrientation', 'landscape');
    hold on;
    subplot(2,2,1);
    hist(nrPossibleEnds, 20);
    xlabel('# possible ends');
    ylabel('occurences');
    title('Histogram of number of possible ends per mission');
    subplot(2,2,2);
    hist([probabilities{:}], 20);
    xlabel('probability');
    ylabel('occurences');
    title('Histogram of probabilities over all missions and possible ends');
    subplot(2,2,[3 4]);
    boxplot(nrObjects, cellfun(@num2str,mat2cell(t,1,ones(1,15)), 'UniformOutput', 0))
    title('Effect of probability thresholding');
    xlabel('proability threshold');
    ylabel('possible ends');
    saveas(gcf, ['C:\Users\mberning\Desktop\problemInspector\' id v(j).name '000Probability.pdf']);
    close all;
end
    
end

