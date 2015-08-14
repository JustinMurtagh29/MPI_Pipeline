function p = changeFileLocationOfParameter(p, newSaveLocation);

    oldSaveLocation = p.saveFolder;
    % Special cases, change data pointers of raw data and currently used CNN
    p.raw.root = strrep(p.raw.root, '/zdata/manuel', '/gaba/u/mberning');
    p.cnn.first = strrep(p.cnn.first, '/zdata/manuel', '/gaba/u/mberning');
    % Recurse over all other elements and replace oldSaveLocation with newSaveLocation
    p = recursiveReplaceInParameter(p, oldSaveLocation, newSaveLocation);

end

function p = recursiveReplaceInParameter(p, oSL, nSL)
    % Only kind of recursive, lots of iteration below xD
    % Get fieldnames of struct
    % Iterate over struct in case it is non scalar
    for i=1:size(p,1)
        for j=1:size(p,2)
            for k=1:size(p,3)
                if ~strcmp(class(p), 'function_handle')
                    switch class(p(i,j,k))
                        case 'struct'
                            fP = fieldnames(p);
                            for f=1:length(fP)
                                p(i,j,k).(fP{f}) = recursiveReplaceInParameter(p(i,j,k).(fP{f}), oSL, nSL);
                            end
                        case 'char'
                            p(i,j,k) = strrep(p(i,j,k), oSL, nSL);
                        case 'cell'
                            for c=1:length(p{i,j,k})
                                p{i,j,k}(c) = recursiveReplaceInParameter(p{i,j,k}(c), oSL, nSL);
                            end
                    end
                end
            end
        end
    end
end

