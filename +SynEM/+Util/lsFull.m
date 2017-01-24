function [ list ] = lsFull( path, args )
%LSFULL Like ls but returns the full path for each file.
% INPUT path: A valid ls input. If only a folder is specified then path
%           should end with a filesep.
%           e.g. /path/to/nmlcollestion/
%                /path/to/data/*.nml
%       args: (Optional on UNIX) String arguments passed to ls.
%           (Default: '');
% OUTPUT list: The same as the ls output but containing the full path to
%           each file.
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

if ispc
    list = ls(path);
    path = [fileparts(path) filesep];
    list = cat(2,repmat(path,size(list,1),1),list);
else
    if exist('args','var') && ~isempty(args)
        args = ['-1', args];
    else
        args = '-1';
    end
    list = ls(path, args);
    list = strsplit(list,'\n');
    list(cellfun(@isempty,list)) = [];
    list = list';
end

end

