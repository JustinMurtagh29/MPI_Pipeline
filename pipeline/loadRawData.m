function raw = loadRawData(root, prefix, bbox, normalize)

raw = readKnossosRoi(root, prefix, bbox);
% Normalize data
if normalize
	raw = normalizeStack(single(raw));
else
	raw = single(raw);
end


end 
