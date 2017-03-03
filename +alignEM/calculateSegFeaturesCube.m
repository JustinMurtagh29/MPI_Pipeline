function calculateSegFeaturesCube(p, cubeNo, fm)
% Calculate SynEM features for SegEM output on a local cube


% load segmentation
bboxFM = bsxfun(@plus, p.local(cubeNo).bboxSmall,[-fm.border', fm.border']./2);
seg = loadSegDataGlobal(p.seg, bboxFM);

load(p.local(cubeNo).segmentFile, 'segments');
interfaces.surface = {segments.PixelIdxList};
interfaces.subseg = cell(1,0); % needed for length check in fm.calculate
clearvars segments

%calculate features
features = fm.calculate(interfaces, seg);

% Save features
if strcmp(fm.mode, 'direction')
    features = features(1:end/2,:); %only save first direction
end
outputFile = [p.local(cubeNo).saveFolder 'InterfaceClassFeatures.mat'];
Util.save(outputFile, features);

end
