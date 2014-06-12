function class = loadClassData(root, prefix, bbox)

% Load data with right border for cnet
class = readKnossosRoi(root, prefix, bbox, 'single', '', 'raw');

end
