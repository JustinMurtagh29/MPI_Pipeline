function job = splitBordersMeta(p)
    mkdir([p.saveFolder 'splitBorder/'])
    functionH = @connectEM.splitBorders;
    inputCell = cellfun(@(x){x}, num2cell(1:50), 'uni', 0);
    cluster = Cluster.getCluster( ...
        '-pe openmp 1', ...
        '-p 0', ...
        '-l h_vmem=64G', ...
        '-l s_rt=23:50:00', ...
        '-l h_rt=24:00:00');
    job = Cluster.startJob( functionH, inputCell, ...
        'sharedInputs', {p},  'sharedInputsLocation', 2, ...
        'name', 'splitBorders', ...
        'cluster', cluster);
    superagglosBorderSplit = {};
    axonLength = {};
    for idx = 1:50
        idx
        temp = load([p.saveFolder 'splitBorder/splitBorder_' num2str(idx)])
        
        superagglosBorderSplit = [superagglosBorderSplit, temp.superagglosBorderSplit(idx:50:end)];
        axonLength = [axonLength, temp.axonLength(idx:50:end)];
    end
    superagglosBorderSplit = cell2mat(superagglosBorderSplit);
    axonLength=cell2mat(axonLength);
    axons = superagglosBorderSplit;
    indBigAxons = axonLength > 5000;
    save(fullfile(p.saveFolder,'aggloState/axons_04.mat'), 'axons','indBigAxons');
    