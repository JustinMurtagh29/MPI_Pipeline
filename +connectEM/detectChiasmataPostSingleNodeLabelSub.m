function detectChiasmataPostSingleNodeLabelSub(startidx, outputFolder)
load([outputFolder 'prep_singlenodelabel.mat']);
for i=startidx:500:length(cc)
    [~, pos{i}, dir{i}, queryIdx{i}] = connectEM.detectChiasmataNodes( ...
        nodes,edges,prob,p,cc{i}(centerOfCC(i))); %#ok
end
save([outputFolder, 'temp_singlenodelabel_' num2str(startidx)],'pos', 'dir', 'queryIdx', 'centerOfCC', 'cc','-v7.3');
end