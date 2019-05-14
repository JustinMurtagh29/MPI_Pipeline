function data = writeJson( file, data )
    file_id = fopen(file, 'w');
    rawData = jsonencode(data);%tojson(data);
    fprintf(file_id, '%s', rawData); % Check github or redo compilation yourself with int64 as integer style
    fclose(file_id);
end

