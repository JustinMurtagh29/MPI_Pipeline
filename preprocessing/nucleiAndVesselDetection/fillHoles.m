function binIm = fillHoles(binIm, mask,holeLoc)
% Dataset specific changes needed here, this function takes mask and
% drill holes into 'outer hull' to be able to use imfill, change drill
% location according to dataset :)
binIm = or(binIm,mask);
if ~exist('holeLoc','var') || isempty(holeLoc)
    holeLoc = round(size(binIm,1)/2);
end
idx = 1;
while binIm(holeLoc,idx) == 1
    binIm(holeLoc,idx) = 0;
    idx = idx + 1;
end

binIm = imfill(binIm, 'holes') & ~mask;
end