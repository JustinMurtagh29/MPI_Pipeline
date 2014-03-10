function output = readKnossosNmlFolder( path )

    % Create a list of all filenames in the folder
    output = struct( 'trace', {}, 'name', {} );
    list = dir( [path '*.nml'] );
    
    % Read every tracing separately and save it in output
    for i = 1 : length( list )
        output(i).name = list(i).name;
        list(i).name = fullfile( path, list(i).name);
        output(i).trace = readKnossosNml( list(i).name, 1 );
    end
end
