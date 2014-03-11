% Startup File manuelCode repo
global GLOBAL_CODE_DIR GLOBAL_DATA_DIR GLOBAL_RESULT_DIR GLOBAL_OUTPUT_DIR GLOBAL_CPU_JM GLOBAL_GPU_JM;
[~,hostname]= system('hostname');
hostnameShort = hostname(1:6);
switch hostnameShort
	case 'fermat'
		base = '/zdata/manuel/';
		GLOBAL_CPU_JM = 'fermat-cpu';
		GLOBAL_GPU_JM = 'fermat-cnn';
		GLOBAL_CODE_DIR = [base 'code' filesep];
		GLOBAL_DATA_DIR = [base 'data' filesep];
	        GLOBAL_RESULT_DIR = [base 'results' filesep];
	        GLOBAL_OUTPUT_DIR = [base 'sync' filesep];
		addpath('/usr/local/jacket/');
		addpath('/usr/local/jacket/engine');
	case 'P1-380'
	        base = 'I:\CortexConnectomics\Manuel\';
		GLOBAL_CPU_JM = 'p1-380-cpu';
		GLOBAL_GPU_JM = 'p1-380-gpu';
		GLOBAL_CODE_DIR = 'C:\manuelCode\';
		GLOBAL_DATA_DIR = [base 'data' filesep];
	        GLOBAL_RESULT_DIR = [base 'results' filesep];
	        GLOBAL_OUTPUT_DIR = [base 'sync' filesep];
	otherwise
		warning(manuelCode:unknownHost,['Warning: Computer ' hostname ' unknown. Please add to startup.m']);
end
clear hostname* base;

% Set path
addpath(genpathGit(GLOBAL_CODE_DIR));

