%% A script to test the close holes algorithm.
load(fullfile('/gaba/u/mberning/results/pipeline/20170217_ROI', 'allParameter.mat'));
graph = load([p.saveFolder 'graphNew.mat'], 'prob', 'edges', 'borderIdx');
meta = load(fullfile(p.saveFolder, 'segmentMeta.mat'), 'segIds', 'point', 'box');
load('/gaba/u/mberning/results/pipeline/20170217_ROI/soma/NucleiCoordinates.mat','rp');
gb = load([p.saveFolder 'globalBorder.mat'], 'borderCoM', 'borderSize');

load('/gaba/u/mberning/results/pipeline/20170217_ROI/aggloState/somas.mat');
voxelSize = [11.24,11.24,28];

aggloSegIds = somas{56,3}
%soma coordinates:  2283, 6035, 2430
somaID = 73

[ newSegIds, name ] = Soma.closeHoles( p, meta, rp, graph, voxelSize, somaID, aggloSegIds )

%% get mappings
segIds = aggloSegIds;
segIds = unique(segIds);
equivalClasses = cell(1,1);
segIds = cat(1, [0], segIds);
equivalClasses{1} = segIds;
script1 = WK.makeMappingScript(max(meta.segIds), equivalClasses, false);

segIds = newSegIds;
segIds = unique(segIds);
equivalClasses = cell(1,1);
segIds = cat(1, [0], segIds);
equivalClasses{1} = segIds;
script2 = WK.makeMappingScript(max(meta.segIds), equivalClasses, false);


%% evaluate
tracing = readNml('/gaba/u/mberning/results/pipeline/20170217_ROI/soma/SomaLocalSurfaceMM_4.nml');

outsideIds = Seg.Global.getSegIds(p, tracing{1}.nodes(:,[1,2,3]));
insideIds = Seg.Global.getSegIds(p, tracing{2}.nodes(:,[1,2,3]));
outsideIds = outsideIds(outsideIds~=0);
insideIds = insideIds(insideIds~=0);

%get correctly added seg ids
notAgglomerated = setdiff(insideIds,aggloSegIds);
correctlyAdded = setdiff(notAgglomerated,newSegIds);
%get wrongly added seg ids
allAdded = setdiff(newSegIds, aggloSegIds);
wronglyAdded = intersect(allAdded, outsideIds);

size(correctlyAdded)
size(wronglyAdded)
