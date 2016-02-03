function subsampleAgain(iF, oF)

    load(iF);
    for i=1:length(issfs)
       issfs{i} = reducepatch(issfs{i}, .1);
    end
    exportSurfaceToAmira(issfs', oF);

end
