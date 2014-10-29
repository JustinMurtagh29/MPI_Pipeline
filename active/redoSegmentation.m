function redoSegmentation( parameter )

aff = loadClassData(parameter.class.root, parameter.class.prefix, parameter.bboxBig);
aff = imcomplement(aff);
map = imimposemin(aff, imextendedmin(aff, .2));
seg = watershed(map, 26);
seg = uint16(seg);

% Save result to KNOSSOS folder
writeKnossosRoi(parameter.seg.root, parameter.seg.prefix, parameter.bboxBig(:,1)', seg, 'uint16');


end

