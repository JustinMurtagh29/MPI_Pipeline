

load('/tmpscratch/mbeining/dendritesSuperPreQuery.mat')                                                            
load('/gaba/u/mbeining/code/Marcel/aboveFiveIndices_superagglos.mat','aboveFiveDend');              

bigDends = superagglos(aboveFiveDend);
smallDends = superagglos(~aboveFiveDend);

connectEM.addQueries(bigDends,'/tmpscratch/mbeining/dendrites2',{'/tmpscratch/kboerg/MBKMB_L4_dendrite_queries_2017_c_nmls/','/tmpscratch/kboerg/MBKMB_L4_dendrite_queries_2017_c_inverse_nmls/'},[],0)

bigDendritesAfterQuery = load('/tmpscratch/mbeining/dendrites2/superagglos.mat');


dendritesPostQuery = [bigDendAfterQuery.superagglos';smallDends];
aboveFiveDend = [true(numel(bigDendAfterQuery.superagglos),1);false(numel(smallDends),1)];
save('/tmpscratch/mbeining/dendritesPostQuery.mat','dendritesPostQuery','aboveFiveDend');