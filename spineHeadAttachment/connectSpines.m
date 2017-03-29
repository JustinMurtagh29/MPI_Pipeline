function connectSpines(blockOf500)
	x = 1;
	save(sprintf('/gaba/u/kboerg/output_spines/touch_%u', blockOf500),'x');
	load('/gaba/u/kboerg/spineHeads_R1_P38.mat');
	load('/gaba/u/kboerg/dendriteIDs.mat');
	load('/gaba/u/kboerg/connM.mat');
	starting = (blockOf500 - 1) * 500 + 1 : min(blockOf500 * 500, length(spineHeads.segIds));
	good_list = connectSpinesSub(spineHeads.segIds(starting), dendriteIDs, connM)
	save(sprintf('/gaba/u/kboerg/output_spines/good_list_%u', blockOf500),'good_list');
end