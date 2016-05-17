function [p, d, e, vE] = determineQuerryLocation(graph, com, cc)
    % Determine querry position, direction and edge based on PCA

    theseCoM = com(cc,:);
    [coeff,score,latent] = princomp(theseCoM);
    [~,minIdx] = min(score(:,1));
    [~,maxIdx] = max(score(:,1));
    p(1,:) = theseCoM(minIdx,:);
    p(2,:) = theseCoM(maxIdx,:);
    d(1,:) = -coeff(:,1);
    d(2,:) = coeff(:,1);
    e(1,:) = determineEdge(graph, cc, minIdx);
    e(2,:) = determineEdge(graph, cc, maxIdx);
    vE = latent./sum(latent);

end

function edge = determineEdge(graph, cc, ccIdx)
    theseN = graph.neighbours{cc(ccIdx)};
    theseP = graph.neighProb{cc(ccIdx)};
    keep = ~ismember(theseN, cc);
    theseN = theseN(keep);
    theseP = theseP(keep);
    if isempty(theseP)
        % Needs to be fixed at some point
        edge = [cc(ccIdx) 0];
    else
        [~,idx] = max(theseP);
        edge = [cc(ccIdx) theseN(idx)];
    end
end


