function target = makeIlastikStack( raw, target )
%MAKEILASTIKSTACK Prepare stack for immport in ilastik.
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

%usual workflow
% load(stacks(1).stackFile)
% target = uint8(target == -1);
% target = flip(rot90(target,3),2);
% makeIlastikStack(raw, target);

[X, Y] = Codat.CNN.Misc.tileTrainingCubes(raw, target, size(target), [0, 0, 0]);
raw = X{1};
target = Y{1};

%write raw
imwrite(raw(:,:,1),'stack.tif');
for z = 2:size(raw,3)
    imwrite(raw(:,:,z),'stack.tif','WriteMode','append');
end

%workaround to get labels in ilastik:
%labels need to be added manually in the .ilp project file (define the
%label in ilastik, trace points at the very end of each parameter range to
%set correct labelset size in ilastik, replace labels with target in .ilp
%file)
% use the commands:
% target = flip(rot90(target,3),2);
% h5disp('data\ilastikProjects\cortexTrainingData.ilp') %display content of
% project file and search for correct stack, then
% h5write('data\ilastikProjects\cortexTrainingData.ilp', '/PixelClassification/LabelSets/labels000/block0000', permute(target,[4 1 2 3]))

% to get the annotated targets from the ilp file:
% targets = h5read('data\ilastikProjects\cortexTrainingData.ilp','/PixelClassification/LabelSets/labels000/block0000');
% targets = flip(rot90(squeeze(targets),3),2);

% to convert the MT labels into the pred array used in movieMakerMMV use
% pred = zeros(100,100,100,3,'single');
% pred(:,:,:,1) = targets == 1;
% pred(:,:,:,2) = targets == 2;
% pred(:,:,:,3) = targets == 3;
end

