function jobMonitor(job)

finished = false;
while ~finished
	clc;
	info = get(job.tasks);
	for i=1:length(info)
		display(['-- Task Nr.: ' num2str(i,'%.2i') ' --']);
		display([info(i).CommandWindowOutput]);
	end
	pause(60);
	states = {info.State};
	isFinished = strcmp('finished', states);
	finished = all(isFinished)
end

end
