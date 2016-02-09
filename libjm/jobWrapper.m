function jobWrapper(fH, iC)

    for i=1:length(iC)
        fH(iC{i}{:});
    end

end

