function mito = mitoDetection( raw, seg )
% Mito detection heuristic for current classifier
% WARNING: segmentation WITH borders needed

% Prepare data
seg = seg ~= 0;
seg = imerode(seg, ones(3, 3, 3));
raw(raw < 50) = 50;
raw(raw > 180) = 180;
raw = (raw - 50) ./ 130;

% Extract Mitos (dark large intracellular regions)
diff = seg - raw;
diff1 = diff > .68;
diff2 = bwareaopen(diff1, 15^3, 26);
mito = imclose(diff2, ones(3,3,3));

end

