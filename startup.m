% Startup File manuelCode repo
global GLOBAL_HOST GLOBAL_CODE_DIR GLOBAL_DATA_DIR GLOBAL_RESULT_DIR GLOBAL_OUTPUT_DIR CLUSTER_CPU CLUSTER_GPU;
[~,hostname]= system('hostname');
% Somehow Matlab gets hostname with weird character at the end
GLOBAL_HOST = hostname(1:end-1);
switch GLOBAL_HOST
    case 'gaba'
        base = '/gaba/u/mberning/';
        GLOBAL_CODE_DIR = [base 'code' filesep];
        GLOBAL_DATA_DIR = [base 'data' filesep];
        GLOBAL_RESULT_DIR = [base 'results' filesep];
        GLOBAL_OUTPUT_DIR = [base 'sync' filesep];
    case 'gaba02'
        base = '/gaba/u/mberning/';
        GLOBAL_CODE_DIR = [base 'code' filesep];
        GLOBAL_DATA_DIR = [base 'data' filesep];
        GLOBAL_RESULT_DIR = [base 'results' filesep];
        GLOBAL_OUTPUT_DIR = [base 'sync' filesep];
    case 'M-01273'
        base = '/home/mberning/localStorage/';
        GLOBAL_CODE_DIR = [base 'code' filesep];
        GLOBAL_DATA_DIR = [base 'data' filesep];
        GLOBAL_RESULT_DIR = [base 'results' filesep];
        GLOBAL_OUTPUT_DIR = [base 'sync' filesep];
    otherwise
        warning(manuelCode:unknownHost,['Warning: Computer ' hostname ' unknown. Please add to startup.m']);
end
% Set path
addpath(genpathGit(GLOBAL_CODE_DIR));
% Log commands
diaryStart = datestr(clock, 30);
diaryName = input('Please enter name for the diary: ', 's');
diary([base 'log' filesep diaryStart '_' diaryName '.log']);
% Set up cluster objects
CLUSTER_CPU = createParcluster('cpu');
CLUSTER_GPU = createParcluster('k40');
% Clear variables not needed anymore
clear hostname base diaryName diaryStart;

