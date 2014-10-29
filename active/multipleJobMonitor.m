function multipleJobMonitor(jobs)

finished = false;
tic;
while ~finished
	clc;
	for i=1:length(jobs)
		info = get(jobs{i}.tasks);
		states = {info.State};
		isFinished = strcmp('finished', states);
		jobFinished(i) = all(isFinished);
	end
	t = toc;
	display(['Time spent: ' num2str(t/3600, '%.2f') ' hours, Jobs finished: ' num2str(sum(jobFinished)) '/' num2str(length(jobFinished))]);
	finished = all(jobFinished);
end

end
