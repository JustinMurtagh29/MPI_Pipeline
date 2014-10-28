% Startup File manuelCode repo
global GLOBAL_HOST GLOBAL_CODE_DIR GLOBAL_DATA_DIR GLOBAL_RESULT_DIR GLOBAL_OUTPUT_DIR GLOBAL_CPU_JM GLOBAL_GPU_JM;
[~,hostname]= system('hostname');
GLOBAL_HOST = hostname(1:end-1);
hostnameShort = hostname(1:4);
switch hostnameShort
	case 'turi'
		base = '/zdata/manuel/';
		GLOBAL_CPU_JM = 'fermat-job-manager';
		GLOBAL_GPU_JM = 'fermat-gpu-manager';
		GLOBAL_CODE_DIR = [base 'code' filesep];
		GLOBAL_DATA_DIR = [base 'data' filesep];
	    GLOBAL_RESULT_DIR = [base 'results' filesep];
	    GLOBAL_OUTPUT_DIR = [base 'sync' filesep];
	case 'ferm'
		base = '/zdata/manuel/';
        GLOBAL_CPU_JM = 'fermat-job-manager';
		GLOBAL_GPU_JM = 'fermat-gpu-manager';
		GLOBAL_CODE_DIR = [base 'code' filesep];
		GLOBAL_DATA_DIR = [base 'data' filesep];
	    GLOBAL_RESULT_DIR = [base 'results' filesep];
	    GLOBAL_OUTPUT_DIR = [base 'sync' filesep];
		addpath('/usr/local/jacket/');
		addpath('/usr/local/jacket/engine');
	case 'P1-3'
	    base = 'I:\CortexConnectomics\Manuel\';
		GLOBAL_CPU_JM = 'p1-380-cpu';
		GLOBAL_GPU_JM = 'p1-380-gpu';
		GLOBAL_CODE_DIR = 'P:\manuelCode\';
		GLOBAL_DATA_DIR = [base 'data' filesep];
	    GLOBAL_RESULT_DIR = [base 'results' filesep];
	    GLOBAL_OUTPUT_DIR = [base 'sync' filesep];
    case 'gaba'
        base = '/gaba/u/mberning/';
		GLOBAL_CPU_JM = 'gabaCPU';
		GLOBAL_GPU_JM = 'gabaGPU';
		GLOBAL_CODE_DIR = [base 'code' filesep];
		GLOBAL_DATA_DIR = [base 'data' filesep];
	    GLOBAL_RESULT_DIR = [base 'results' filesep];
	    GLOBAL_OUTPUT_DIR = [base 'sync' filesep];
 	case 'quak'
	  base = '/zdata/manuel/';
    GLOBAL_CPU_JM = 'fermat-job-manager';
		GLOBAL_GPU_JM = 'fermat-gpu-manager';
		GLOBAL_CODE_DIR = [base 'code' filesep];
		GLOBAL_DATA_DIR = [base 'data' filesep];
	  GLOBAL_RESULT_DIR = [base 'results' filesep];
	  GLOBAL_OUTPUT_DIR = [base 'sync' filesep];
  otherwise
		warning(manuelCode:unknownHost,['Warning: Computer ' hostname ' unknown. Please add to startup.m']);
end
clear hostname* base;

% Set path
addpath(genpathGit(GLOBAL_CODE_DIR));
% Load parameters for big pipeline test
[p,pT] = setParameterSettingsBig('20140312T141921');
% Log commands
diaryStart = datestr(clock, 30);
diaryName = input('Please enter name for the diary: ', 's');
diary(['/zdata/manuel/log/' diaryStart diaryName '.log']);

