function predictLocalCube(saveFolder, classifierFile, outputFilenameInLocalFolder)
    if ~exist('classifierFile', 'var') || isempty(classifierFile)
        classifierFile = '/gaba/u/sahilloo/Mk1_F6_JS_SubI_v1/connect-em/edgeClassifier/20190329T183128.mat';
    end
    
    if ~exist('outputFilenameInLocalFolder', 'var') ...
            || isempty(outputFilenameInLocalFolder)
        outputFilenameInLocalFolder = 'neuriteContinuityProb.mat';
    end
    
    % Load needed data
    classifier = load(classifierFile);
    classifier = classifier.classifier;
    
    borders = load(fullfile(saveFolder, 'borders.mat'));
    borders = borders.borders;
    
    rawF = load(fullfile(saveFolder, 'InterfaceRawFeatures.mat'));
    classF = load(fullfile(saveFolder, 'InterfaceClassFeatures.mat'));
    features = cat(2, rawF.features, classF.features);
    clear rawF classF;
    
    % set all edges that are not classified to 0 probability
    classifiedBorderIdx = cat(1, borders.Area) > 10;
    prob = zeros(size(classifiedBorderIdx));
    
    if ~isempty(features)
        [~, scores] = classifier.predict(features);
        sigmoid = @(x) 1 ./ (1 + exp(-x));
        prob(classifiedBorderIdx) = sigmoid(scores(:, 1));
    end
    
    outFile = fullfile(saveFolder, outputFilenameInLocalFolder);
    Util.save(outFile, prob);
end

