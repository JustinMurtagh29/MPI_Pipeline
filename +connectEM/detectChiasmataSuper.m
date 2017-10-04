function detectChiasmataSuper(startidx,p, useSphereClustering)

    if nargin < 3
        useSphereClustering = false;
    end
    temp = load(fullfile(p.saveFolder,'aggloState/axons_06_c.mat'));
    temp.axons = temp.axons(temp.indBigAxons);
    numstr='36';
    for idx = startidx : 500 : length(temp.axons)
        outputFolder = fullfile('/tmpscratch/kboerg/chiasmata', ['chiasmataX', numstr, '_' num2str(floor(idx/100)) '/visX', numstr, '_' num2str(idx) '/']);
        mkdir(outputFolder);
        if useSphereClustering
            connectEM.detectChiasmataSphereClustering(p, temp.axons(idx).nodes(:, 1:3), temp.axons(idx).edges, false, outputFolder);
        else
            connectEM.detectChiasmata(p, temp.axons(idx).nodes(:, 1:3), temp.axons(idx).edges, false, outputFolder);
        end
    end
end
