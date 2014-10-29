function miniBorder(parameter)

jm = findResource('scheduler', 'type', 'jobmanager', 'configuration', 'fermat-cpu');
job = createJob(jm, 'configuration', 'fermat-cpu');
load(parameter.edgeFile);

i = 0;
while i < size(edges,1)
	batch = (i+1):min(i+1000,size(edges,1));
	inputCell = {parameter, edges(batch,:)};
	createTask(job, @borderCalculation, 3, inputCell, 'configuration', 'fermat-cpu');
	i = batch(end);
end

% Start Job
submit(job)

% Postprocessing
waitForState(job);
data = getAllOutputArguments(job);
edges = cat(1,data{:,1});
border = cat(2,data{:,2});
weights = cat(1,data{:,3});
save(parameter.edgeFile, 'edges');
save(parameter.borderFile, 'border');
save(parameter.weightFile, 'weights');
destroy(job);

end

