function ekSegmentationToJm()

xRange = 1:34; % Segmentiert cubes 2-35 (leaves out boundary cubes & takes into consideration cube number starting at 0)
yRange = 1:41; % Segmentriert cubes 2-42
dir = '/zdata/manuel/results/CNNfwdPass/14-Feb-2013618765/';
jm = findResource('scheduler', 'configuration', 'local_1');
job = createJob(jm, 'configuration', 'local_1');
for k=38:44
	inputargs = {dir, xRange,yRange,k-1};
	createTask(job, @segment_ek_0563ForPaper, 0, inputargs, 'configuration', 'local_1');
end
submit(job);

end

