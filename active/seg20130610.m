function seg20130610( parameter, bbox )

cube = floor((bbox(:,1)-1)/128); 
bbox = bbox + [-10 10; -10 10; -10 10];

aff = loadClassData(parameter.class.root, parameter.class.prefix, bbox);
aff = imcomplement(aff);
fgm1 = imextendedmin(aff, .39, 26);
fgm1 = bwareaopen(fgm1, 10);
map = imimposemin(aff, fgm1);
seg = watershed(map, 26);
seg = uint16(seg);

% Save result to KNOSSOS folder
writeKnossosCube(parameter.seg.root, parameter.seg.prefix, cube, seg, 'uint16');


end

