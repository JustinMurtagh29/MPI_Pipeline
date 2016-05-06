function trainingData = getTrainingData( )
%GETTRAININGDATA Load training data from standard location.
% OUTPUT trainingData: struct containing the training data


global GLOBAL_DATA_DIR
trainingFolderPath = [GLOBAL_DATA_DIR 'Features' filesep 'Training' filesep];
testFolderPath = [GLOBAL_DATA_DIR 'Features' filesep 'Test' filesep];
[x_train,y_train,featureMap,trainCubeNames] = getTrainingDataFrom(trainingFolderPath);
[x_test,y_test,~,testCubeNames,intLookupTest] = getTrainingDataFrom(testFolderPath);
trainingData.x_test = x_test;
trainingData.y_test = y_test;
trainingData.x_train = x_train;
trainingData.y_train = y_train;
trainingData.featureMap = featureMap;
trainingData.testCubeNames = testCubeNames;
trainingData.trainCubeNames = trainCubeNames;
trainingData.lookupTest = intLookupTest;

end

