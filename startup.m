% Startup File manuelCode repo
global GLOBAL_CODE_DIR GLOBAL_DATA_DIR GLOBAL_RESULT_DIR GLOBAL_OUTPUT_DIR;
[~,hostname]= system('hostname');
hostnameShort = hostname(1:6);
switch hostnameShort
	case 'fermat'
		base = '/zdata/manuel/';
		GLOBAL_CODE_DIR = [base 'code' filesep];
		GLOBAL_DATA_DIR = [base 'data' filesep];
	        GLOBAL_RESULT_DIR = [base 'results' filesep];
	        GLOBAL_OUTPUT_DIR = [base 'sync' filesep];
	case 'P1-380'
	        base = 'I:\CortexConnectomics\Manuel\';
		GLOBAL_CODE_DIR = 'C:\manuelCode\';
		GLOBAL_DATA_DIR = [base 'data' filesep];
	        GLOBAL_RESULT_DIR = [base 'results' filesep];
	        GLOBAL_OUTPUT_DIR = [base 'sync' filesep];
	otherwise
		warning(manuelCode:unknownHost,['Warning: Computer ' hostname ' unknown. Please add to startup.m']);
end
clear hostname* base;

