function job = makeInterfacePRedictionsOnCluster(p)

display('Randomly sampling 100 cubes to calulate feature quantiles for Interface Classification...');
calculateFeatureQuantilesForIC(p);
display('Finished calculating feature quantiles! Now making predictions');

functionH = @makeInterfacePredictions;
inputCell = {p};

job=startCPU(functionH,inputCell,'InterfacePredictions');

end
