function detectChiasmataSuper(startidx,p)
    temp = load(fullfile(p.saveFolder,'aggloState/axons_0X.mat'));
    for idx = startidx : 500 : length(temp.axons)
        outputFolder = fullfile(p.saveFolder, 'chiasmata' ['chiasmataX', numstr, '_' num2str(floor(idx_agglo/100)) '/']);
        mkdir(outputFolder);
        detectChiasmata(p, temp.axons(1).nodes, temp.axons(1).edges, false, outputFolder);
    end
end