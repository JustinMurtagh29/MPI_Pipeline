function result = runCNNEvaluation( saveFolder, cnnFile )
%RUNCNNEVALUATION Run the whole segmentation evaluation pipeline.
% INPUT saveFolder: Root folder where all segmentations are stores
%           (see setParameters07x2)
%       cnnFile: Full path to cnn.
% OUTPUT result: see parameterSearchSeg
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

%create segmenation folder
[~,cnnName,~] = fileparts(cnnFile);
[p,pT] = Codat.Segmentation.setParameters07x2(saveFolder, cnnName, cnnName);

%save cnn to folder
copyfile(cnnFile,[p.cnn.saveFile, '.mat']);

%run cnn
fprintf('[%s] runCNNEvaluation - Starting membrane prediction.\n',...
    datestr(now));
job = Codat.Segmentation.minicubeFwdPass(pT);

%wait for output
fprintf('[%s] runCNNEvaluation - Waiting for membrane prediction.\n', ...
    datestr(now));
wait(job);

%run parameter search
fprintf('[%s] runCNNEvaluation - Running segmentation parameter search.\n', ...
    datestr(now));
result = Codat.Segmentation.parameterSearchSeg(pT);

%save result in segmentation folder
filename = [pT.saveFolder, 'result.mat'];
fprintf('[%s] runCNNEvaluation - Saving result to %s.\n', ...
    datestr(now), filename);
save(filename, 'result');

end

