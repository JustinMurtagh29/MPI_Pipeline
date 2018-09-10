
%Alessandro magic


synSegments = cellfun(@vertcat, syn.synapses.presynId, syn.synapses.postsynId, 'UniformOutput', false);

%maxSegId = Seg.Global.getMaxSegId(param);
maxSegId = 15030572;

lut = false(maxSegId, 1);
lut(setdiff(seg, 0)) = true;

tic
synInBox = cellfun(@(synSegIds) any(lut(synSegIds) == true), syn.synapses.postsynId);
toc


% synInBox = false();
% for i = 1:numel(synSegments)
%     synInBox(i) = any(lut(synSegment{i}));
% end