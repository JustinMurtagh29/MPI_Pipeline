function detectChiasmataSuperAxons_05b(startidx,p)
    temp = load(fullfile(p.saveFolder,'aggloState/axons_05_b.mat'));
    temp.axons = temp.axons(temp.indBigAxons);
    numstr='32';
    for idx = startidx : 500 : length(temp.axons)
        outputFolder = fullfile(p.saveFolder, 'aggloState/chiasmata/chiasmata_axons_05b', ['chiasmataX', numstr, '_' num2str(floor(idx/100)) '/visX', numstr, '_' num2str(idx) '/']);
        mkdir(outputFolder);
        connectEM.detectChiasmata(p, temp.axons(idx).nodes(:, 1:3), temp.axons(idx).edges, false, outputFolder);
    end
end
