function detectChiasmataPostSingleNodeLabelSub(startidx)
load(['/tmpscratch/kboerg/visX19_0/visX19_1/prep_singlenodelabel.mat']);
for i=startidx:500:length(cc)
    [pos, dir, queryIdx] = connectEM.detectChiasmataPostSingleNodeLabelSubSub(nodes,edges,prob,p,cc,centerOfCC,pos,dir,queryIdx,i)
end
save(['/tmpscratch/kboerg/visX19_0/visX19_1/temp_singlenodelabel_' num2str(startidx)],'pos', 'dir', 'queryIdx', 'centerOfCC', 'cc','-v7.3');
end
function dummy()
end