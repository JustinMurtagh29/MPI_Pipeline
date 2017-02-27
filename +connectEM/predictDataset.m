function job = predictDataset(p);

functionHandle = @connectEM.predictLocalCube;
saveFolders = {p.local(:).saveFolder};
saveFolders = cellfun(@(x){x}, saveFolders, 'uni', 0);
job = startCPU(functionHandle, saveFolders, 'neuriteContinuityPrediction');

end

