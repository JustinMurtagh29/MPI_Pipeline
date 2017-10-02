function detectChiasmataSuper(startidx,p)
    temp = load(fullfile(p.saveFolder,'aggloState/axons_06_c.mat'));
    temp.axons = temp.axons(temp.indBigAxons);
    numstr='34';
    for idx = startidx : 500 : length(temp.axons)
        outputFolder = fullfile('/tmpscratch/kboerg/chiasmata', ['chiasmataX', numstr, '_' num2str(floor(idx/100)) '/visX', numstr, '_' num2str(idx) '/']);
        mkdir(outputFolder);
        connectEM.detectChiasmata(p, temp.axons(idx).nodes(:, 1:3), temp.axons(idx).edges, false, outputFolder);
    end
end
