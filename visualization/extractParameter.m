function parameter = extractParameter( data )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
parameter = struct([]);
for i=1:size(data,1)
	for j=1:size(data,2)
        if ~isempty(data(i,j).directory)
            parameter(i,j).name = strrep(strrep(data(i,j).directory(1:end-1), 'D:\cnetRetinaForPaper\CNNtrain\', ''), '\', '-');
            load([data(i,j).directory data(i,j).files{1}]);
            parameter(i,j).nrFM = cnet. numFeature;
            parameter(i,j).sizeFM = prod(cnet.filterSize);
            parameter(i,j). nrHiddenLayer = cnet.numHiddenLayer;
            parameter(i,j).nrOutputs = cnet.numLabels;
            if ~isempty(cnet.run.etaW)
                parameter(i,j).etaW = cnet.run.etaW(1);
                parameter(i,j).etaB = cnet.run.etaB(1);
                parameter(i,j).tauW = (cnet.run.etaW(1)-cnet.run.etaW(2))/cnet.run.etaW(1);
                parameter(i,j).tauB = (cnet.run.etaB(1)-cnet.run.etaB(2))/cnet.run.etaB(1);
            else
                parameter(i,j).etaW = cnet.run.wStart;
                parameter(i,j).etaB = cnet.run.bStart;
                parameter(i,j).tauW = (cnet.run.wStart-cnet.run.wEnd)/cnet.run.wStart;
                parameter(i,j).tauB = (cnet.run.wStart-cnet.run.bEnd)/cnet.run.bStart;
            end
            parameter(i,j).batchSize = prod(cnet.outputSize);
            parameter(i,j).iterations = prod(cnet.run.maxIterMini);
            parameter(i,j).directory = data(i,j).directory;
            parameter(i,j).files = data(i,j).files;
        end
	end
end

end

