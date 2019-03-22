function runBox(p, bbox,saveFolder, border, options, cnn17, cnn17_2, cnn17_4})
% make svm predictions
%
% Author: Sahil Loomba <sahil.loomba@brain.mpg.de>

info = Util.runInfo();

%% predictions using vanilla CNNs
raw = Seg.IO.loadRaw(p, bbox, border ./ 2, true);
raw = uint8(raw*22 + 122);
raw = single(raw) ./ 255;

% cnn_17
pred17 = Codat.CNN.Misc.predictCube(cnn17, raw, options);

% cnn17_2
pred17_2 = Codat.CNN.Misc.predictCube(cnn17_2, raw, options);

% cnn17_4
pred17_4 = Codat.CNN.Misc.predictCube(cnn17_4, raw, options);

Util.save(fullfile(saveFolder, 'svmPredictions.mat'),...
             pred17, pred17_2, pred17_4, bbox, info);

end
