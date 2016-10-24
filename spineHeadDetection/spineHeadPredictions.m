
function [spineHeadsPerUMcube,spineHeads]=spineHeadPredictions(p,th)
%To count the number of spine heads detected. Threshold can be set via variable 'th'
%th=0;
spineHeadCounterTotal=0;
spineHeadCounter=zeros(numel(p.local),1);

load([p.saveFolder 'globalCoMList.mat']); % loads globalComList

spineHeads.Ids=[];
spineHeads.CoMs=[];

for idx=1:numel(p.local)

	m=load([p.local(idx).saveFolder 'spineHeadScores.mat']);
	allScores=m.spineHeadScores(:,3);
	allIds = m.spineHeadScores(:,1);
	%Predict spine head ids using threshold 'th'
	predictedIds = allIds(allScores(:)>th);
	predictedCoMs = globalCoMList(predictedIds,:);
	%store the predicted ids and CoMs in global struct
	spineHeads.Ids=cat(1,spineHeads.Ids,predictedIds);
	spineHeads.CoMs=cat(1,spineHeads.CoMs,predictedCoMs);
	%Count numbre of spine heads in local cube
	spineHeadCounter(idx)=length(predictedIds);

end

% To calculate spine head denstiy per um cube
spineHeadCounterTotal=sum(spineHeadCounter(:));

voxels=p.bbox(:,2)-p.bbox(:,1);
voxelSize=p.raw.voxelSize;
volumeNM=prod(voxelSize.*voxels');
volumeUM=volumeNM/1.0000e+09;

spineHeadsPerUMcube=spineHeadCounterTotal/volumeUM;

save([p.saveFolder 'spineHeads.mat'],'spineHeads')

end
