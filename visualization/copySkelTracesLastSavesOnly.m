function copySkelTracesLastSavesOnly( oldPath, newPath )

% COPYSKLETRACESLASTSAVESONLY: Copy only highest number of each file prefix
% from oldPath to newPath
%   
%   The function has the following arguments:
%       OLDPATH: Give the root directory of all tracing files, e.g.
%           'E:\e_k0563\All_ek0563_Skeletons\'. Don't omit the "\"!
%       NEWPATH: Give the target directory, e.g. 
%           'E:\e_k0563\All_ek0563_Skeletons_LastSavesOnly\'. Don't omit the "\"!
%       
%   => copySkelTracesLastSavesOnly( 'E:\e_k0563\All_ek0563_Skeletons\', 'E:\e_k0563\All_ek0563_Skeletons_LastSavesOnly\' )
%

    % Collect all the traces that are present twice (once with -1 or _1 suffix)
    listDouble = dir( [oldPath '*-1.nml'] );
    listDouble = [ listDouble; dir( [oldPath '*_1.nml'] ) ];
    
    % Collect all filenames of all .nml files
    list = dir( [oldPath '*.nml'] );

    % Copy the filenames from the listDouble struct into a cell array
    cellArrayDouble = cell( size( listDouble, 1 ), 1 );
    
    for i = 1 : size( listDouble, 1 )
        cellArrayDouble{i} = listDouble(i).name;
    end

    % Copy all the filenames (without numbers) from the list struct into a cell array, but
    % leave out the redundant files (those filenames saved in cellArrayDouble)
    cellArrayOfNames = cell( size( list, 1 ), 1 );
    
    for i = 1 : size( list, 1 )
        if ~(strcmp(list(i).name, cellArrayDouble))
            cellArrayOfNames{i} = list(i).name( 1 : end - 8 );
        end
    end
    
    % Remove empty cells from cellArrayOfNames
    cellArrayOfNames = cellArrayOfNames( ~cellfun( 'isempty', cellArrayOfNames ) );
    
    % Find the unique file names in cellArrayOfNames
    tracings = unique( cellArrayOfNames );

    % Copy files from the destination folder to the target
    for i = 1 : size( tracings, 1 )
        
        % Find all the files with the unique filename saved in tracings
        listThisTracing = dir( [oldPath tracings{i} '*.nml'] );
        cellArrayTracing = cell( size( listThisTracing, 1 ), 1 );
        
        % Save only the numbers of this tracing into cellArrayTracing
        for j = 1 : size( listThisTracing, 1 )
            cellArrayTracing{j} = str2double( listThisTracing(j).name( end - 6 : end - 4 ) );
        end
        
        % Determine the maximum number of this tracing and use it to find
        % the last save to copy to the target
        maxTracing = max( cell2mat( cellArrayTracing ) );
        copyfile( [oldPath tracings{i} '.' num2str( maxTracing, '%03i' ) '.nml'], ...
        [newPath  tracings{i} '.' num2str( maxTracing, '%03i' ) '.nml'] );
    end
end