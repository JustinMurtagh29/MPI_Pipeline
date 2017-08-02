function detectChiasmataPostSingleNodeLabelSub(startidx)
load(['/tmpscratch/kboerg/visX19_0/visX19_1/prep_singlenodelabel.mat']);
for i=startidx:5000:length(cc)
    [pos, dir, queryIdx] = detectChiasmataPostSingleNodeLabelSubSub(node,edges,prob,p,cc,centerOfCC,pos,dir,queryIdx,i)
end
save([outputFolder 'temp_singlenodelabel_' num2str(startidx)],'pos', 'dir', 'queryIdx', 'centerOfCC', 'cc','-v7.3');
end
function dummy()
end