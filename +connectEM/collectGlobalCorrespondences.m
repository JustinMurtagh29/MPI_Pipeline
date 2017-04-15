function corrEdges = collectGlobalCorrespondences(p)

    saveFolder = [p.saveFolder 'correspondencesNew' filesep];
    files = dir([saveFolder '*.mat']);
    
    corrEdges = zeros(0,2,'uint32');
    for i=1:length(files)
        load([saveFolder files(i).name]);
        countsCnormMin = countsC ./ min(arrayfun(@(x)countsS(uniqueSegments == x), uniqueCorrespondences),[],2);
        theseEdges = uniqueCorrespondences(countsC > 10 & countsCnormMin > 0.9,:);
        corrEdges(end+1:end+size(theseEdges,1),:) = theseEdges;
    end
end

