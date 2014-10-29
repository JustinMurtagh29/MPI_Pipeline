function miniEdges(parameter)

jm = findResource('scheduler', 'type', 'jobmanager', 'configuration', 'fermat-cpu');
job = createJob(jm, 'configuration', 'fermat-cpu');
seg = loadSegData(parameter.seg.root, parameter.seg.prefix, parameter.bboxSmall);
ids = unique(seg);
ids(ids == 0) = [];

i = 0;
while i < length(ids)
	batch = (i+1):min(i+1000,length(ids));
	inputCell = {parameter, ids(batch)};
	createTask(job, @findEdges, 1, inputCell, 'configuration', 'fermat-cpu');
	i = batch(end);	
end
for i=1:length(parameter.filter)
	createTask(job, @filter3d, 0, {parameter 'aff' i}, 'configuration', 'fermat-cpu');
end
for i=1:length(parameter.filter)
	createTask(job, @filter3d, 0, {parameter 'raw' i}, 'configuration', 'fermat-cpu');
end

% Start Job
submit(job)

% Postprocessing
waitForState(job);
data = getAllOutputArguments(job);
edges = cat(1,data{1:end-12,1});
edges = sort(cell2mat(edges),2);
edges = unique(edges, 'rows');
save(parameter.edgeFile, 'edges');
destroy(job);

end

