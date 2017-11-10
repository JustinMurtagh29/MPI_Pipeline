function nrExits = detectChiasmataSub(startidx, outputFolder)
    load(fullfile(outputFolder, 'prep.mat'));
    nrExits = zeros(size(nodes, 1), 1);
    for i=startidx:5000:size(nodes,1)
        nrExits(i) = connectEM.detectChiasmataNodes( ...
            nodes, edges, ones(size(edges,1), 1), p, i);
    end
    save([outputFolder 'temp_' num2str(startidx)], 'nrExits', '-v7.3');
end
