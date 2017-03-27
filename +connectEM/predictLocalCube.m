function predictLocalCube(saveFolder, classifierFile, outputFilenameInLocalFolder)

    if nargin < 3
        outputFilenameInLocalFolder = 'neuriteContinuityProb.mat';
    end

    if nargin < 2
        classifierFile = '/gaba/u/mberning/results/edgeClassifier/20170322T153247.mat';
    end
    if exist([saveFolder 'InterfaceRawFeatures.mat'],'file')
       % Load needed data
       load(classifierFile);
       load([saveFolder 'borders.mat']);
       rawF = load([saveFolder 'InterfaceRawFeatures.mat']);
       classF = load([saveFolder 'InterfaceClassFeatures.mat']);
       features = cat(2, rawF.features, classF.features);
       clear rawF classF;
   
       % set all edges that are not classified to 0 probability
       classifiedBorderIdx = cat(1,borders(:).Area) > 10;
       prob = zeros(length(classifiedBorderIdx), 1);
       [~,scores] = classifier.predict(features);
       sigmoid = @(x)1./(1+exp(-1.*x));
       prob(classifiedBorderIdx) = sigmoid(scores(:,1));
       save([saveFolder outputFilenameInLocalFolder], 'prob');
    end
end

