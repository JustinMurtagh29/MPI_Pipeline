function compareChiasmataDetections(p)

    axons = load(fullfile(p.saveFolder,'aggloState/axons_06_c.mat'));
    idxBig = find(axons.indBigAxons)';
    axons = axons.axons(idxBig);

    for idx = 1:length(idxBig)
        idx_agglo = idxBig(idx);
        result1 = load(['/tmpscratch/kboerg/chiasmata/chiasmataX34_' num2str(floor(idx_agglo/100)) '/visX34_' num2str(idx_agglo) '/result.mat']);
        result1 = result1.output;
        result2 = load(['/tmpscratch/kboerg/chiasmata/chiasmataX35_' num2str(floor(idx_agglo/100)) '/visX35_' num2str(idx_agglo) '/result.mat']);
        result2 = result2.output;
        nrChiasmataNodes1(idx) = sum(result1.isIntersection);
        nrChiasmataNodes2(idx) = sum(result2.isIntersection);
        nrChiasmataNodesAdded(idx) = sum(~result1.isIntersection & result2.isIntersection);
    end

    display('Number chiasmata nodes part of sphere hull approach but not graph topology based approach:');
    display(num2str(sum(nrChiasmataNodesAdded)));

end

