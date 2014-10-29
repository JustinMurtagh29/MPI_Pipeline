function miniFwdPass(parameter)

jm = findResource('scheduler', 'type', 'jobmanager', 'configuration', 'fermat-cpu');
job = createJob(jm, 'configuration', 'fermat-cpu');
% Load last CNN from training batch given in parameter
load([parameter.cnn.root 'saveNet0000000001.mat'], 'cnet');
cnet = cnet.loadLastCNN;
cnet.run.actvtClass = @single;
% Start a task for each cube within bbox, on CPU
for x = parameter.cubesBig(1,1):parameter.cubesBig(1,2)
	for y = parameter.cubesBig(2,1):parameter.cubesBig(2,2)
		for z = parameter.cubesBig(3,1):parameter.cubesBig(3,2)
			bbox(:,1) = [x y z].*128+1;
			bbox(:,2) = ([x y z]+1)*128;
			inputCell = {cnet, parameter.raw, parameter.class, bbox};
			createTask(job, parameter.class.func, 0, inputCell, 'configuration', 'fermat-cpu');
		end
	end
end
% Start job
submit(job);
waitForState(job);
destroy(job);

end

