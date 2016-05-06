result = [];
%%
result = run(Codat.NN.Tests.testConvFcts);
result = cat(2,result,run(testCase));
%%
result = run(Codat.NN.Tests.testBatchNormalization);
result = cat(2,result,run(testCase));
%% result
rt = table(result)