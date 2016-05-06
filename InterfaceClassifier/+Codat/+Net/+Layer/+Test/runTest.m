result = [];

%% gradient test
testCase = Codat.Net.Layer.Test.LayerNumericalGradientTest();
result = cat(2,result,run(testCase));

%% convolutional layer test
testCase = Codat.Net.Layer.Test.denseConvLayerTest();
result = cat(2,result,run(testCase));

%% deconvolutional layer test
testCase = Codat.Net.Layer.Test.denseDeconvLayerTest();
result = cat(2,result,run(testCase));

%% result
rt = table(result)