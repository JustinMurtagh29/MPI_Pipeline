function borderCalc(idx)
load('/gaba/u/mberning/results/pipeline/20170217_ROI/allParameterWithSynapses.mat');
load('/gaba/u/mberning/results/pipeline/20170217_ROI/aggloState/dendrites_09.mat');
load('/gaba/u/mberning/results/pipeline/20170217_ROI/aggloState/axons_08_a.mat');
seg = readKnossosRoi('/gaba/u/mberning/results/pipeline/20170217_ROI/globalSeg/','2012-09-28_ex145_07x2_ROI2016_corrected_mag1', p.local(idx).bboxSmall, 'uint32');

filterNan = @(x)x(~isnan(x));

dendriteNodes = cat(1, dendrites(indBigDends).nodes);
dendriteIDs = dendriteNodes(:,4);
dendriteLookup(filterNan(dendriteIDs))=repelem(find(indBigDends),arrayfun(@(x)sum(~isnan(x.nodes(:,4))), dendrites(indBigDends)));

axonNodes = cat(1, axons(indBigAxons).nodes);
axonIDs = axonNodes(:,4);
axonLookup(filterNan(axonIDs))=repelem(find(indBigAxons),arrayfun(@(x)sum(~isnan(x.nodes(:,4))), axons(indBigAxons)));
%axon ids will be negative
axonLookup=-axonLookup;

generalLookup = dendriteLookup;
generalLookup(axonLookup~=0)=axonLookup(axonLookup~=0);
if max(seg(:)) > length(generalLookup)
    generalLookup(max(seg(:))) = 0;
end
seg2 = double(seg);
seg2(seg2~=0)=generalLookup(seg2(seg2~=0));
[edges, ind] = connectEM.borders.codeBenedikt(seg2);
[X,Y,Z] = ind2sub(size(seg2)+2,ind); 
%the +2 and -1 is for fixing the padding
indSmall = sub2ind(ceil(size(seg2)./[32,32,16]),ceil((X-1)/32),ceil((Y-1)/32),ceil((Z-1)/16));
for idx2 = 1 : prod(ceil(size(seg2)./[32,32,16]))
    [C,ia,ic] = unique(edges(indSmall==idx2,:), 'rows');
    freqs = tabulate(ic);
    if ~isempty(freqs)
        findings{idx2} = [C, freqs(:,2)];
    else
        findings{idx2} = [];
    end
end
save(['/tmpscratch/kboerg/borders2/borders_' num2str(idx)],'findings');

