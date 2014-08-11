function [ output_args ] = gliaClusterTest( input_args )


pipeline(p,pT,'graphConstructionLR');

pipeline(p,pT,'graphConstruction');

pipeline(p,pT,'filterbankLR');

pipeline(p,pT,'filterbank');

pipeline(p,pT,'prepareTrainingData Glia');

pipeline(p,pT,'prepareGP Glia');

pipeline(p,pT,'applyGP Glia');


end

