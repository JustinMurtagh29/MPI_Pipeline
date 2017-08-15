function axonQueryAnalysisSub(startidx)
temp = load('/gaba/scratch/mberning/axonQueryGeneration/beforeQueryGeneration.mat', 'axonsNew');
axons = temp.axonsNew;
load('/gaba/scratch/mberning/edgesGTall.mat');

for idx = startidx : 500 : length(axons)
    usededges{idx} = find(all(ismember(edgesGTall,axons{idx}),2));
end
save(['/tmpscratch/kboerg/20170810axonQueryAnalysis/output' num2str(startidx)],'usededges');
