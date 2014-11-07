function data = readJson( file )
    file_id = fopen(file);
    rawData = fscanf(file_id, '%c');
    data = fromjson(rawData); % Check github or redo compilation yourself with int64 as integer style
    fclose(file_id);
end

