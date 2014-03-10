% Startup File Manuel Berning

[~,hostname]= system('hostname');
hostname = hostname(1:6);
switch hostname
	case 'fermat'
		codeDir = '/zdata/manuel/code/';
		dataDir = '/zdata/manuel/data/'; 
		jm = findResource('scheduler', 'type', 'jobmanager', 'LookupURL', 'fermat01');
		cpuJobs	= findJob(jm(1));
		gpuJobs = findJob(jm(2));
		display(jm);
	case 'P1-377'
		codeDir = '/home/mberning/code/';
		dataDir = '/data/';	
	case 'P1-380'
    		codeDir = 'C:\code\';
		dataDir = 'C:\data\';
	otherwise
		warning(['Warning: Computer ' hostname ' unknown. ']);
		return;
end
addpath(genpath(codeDir));
% In order to not shadow matlab functions, take out KLEE
rmpath([codeDir 'KLEE']);

