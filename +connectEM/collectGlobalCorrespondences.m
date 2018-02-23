function [corrEdges, corrCom] = collectGlobalCorrespondences(p)
% Author: Manuel Berning <manuel.berning@brain.mpg.de>
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

    saveFolder = p.correspondence.saveFolder;
    files = dir([saveFolder '*.mat']);
    
    corrEdges = zeros(0,2,'uint32');
    corrCom = zeros(0,3,'uint32');
    for i=1:length(files)
        load([saveFolder files(i).name]);
        countsCnormMin = countsC ./ min(arrayfun(@(x)countsS(uniqueSegments == x), uniqueCorrespondences),[],2);
        theseEdges = uniqueCorrespondences(countsC > 10 & countsCnormMin > 0.9,:);
        if nargout > 1
            theseComs = com(countsC > 10 & countsCnormMin > 0.9, :);
        end
        corrEdges(end+1:end+size(theseEdges,1),:) = theseEdges;
        corrCom(end+1:end+size(theseEdges,1),:) = theseComs;
    end
end

