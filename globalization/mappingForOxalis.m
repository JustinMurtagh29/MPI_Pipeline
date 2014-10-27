function mappingForOxalis(p)

% Implement better at some point
load([p.seg.root 'numEl.mat']);
load([p.seg.root 'globalMapping.mat']);

% Find maximum id
load(p.local(end,end,end).segFile);
nrGlobalIDs = max(seg(:));
numEl = numElTotal(end,end,end) + uint32(nrGlobalIDs);

% Create a vector of zeros with the right length
mapping = zeros(numEl, 1, 'uint32');

% Replace mapped IDs with first ID
for i=1:length(components)
    for j=2:length(components{i})
        mapping(components{i}(j)) = components{i}(1);
    end
end

fid = fopen([p.seg.root '2012-09-28_ex145_07x2_corrected_segmentation_mapping.map'], 'w');
fwrite(fid, mapping, 'uint32');
fclose(fid);

end
