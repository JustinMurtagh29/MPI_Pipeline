function [weights, weightNames] = shapeFeatures(segments,siz)

    for i=1:length(segments)
        % switch from ind to sub
        [x,y,z] = ind2sub(siz,double(segments(i).PixelIdxList'));
        PixelList = [x y z];             
        objSize = length(PixelList);
        
        [COEFF,SCORE,latent] = princomp(PixelList);
        pcFracs = latent./sum(latent);
        
        weights(i,1) = log(objSize);                     
        weights(i,2) = pcFracs(1);
        weights(i,3) = pcFracs(2);
        weights(i,4) = pcFracs(3);
        weights(i,5) = pcFracs(3)/pcFracs(2);
        weights(i,6) = pcFracs(3)/pcFracs(1);
        weights(i,7) = pcFracs(2)/pcFracs(1);
        weights(i,8) = -pcFracs(1)*log(pcFracs(1))-pcFracs(2)*log(pcFracs(2))-pcFracs(3)*log(pcFracs(3)); % entropy
        %work on that.. cant handle small segments, bad calculations...
        %weights(i,9) = calcConcavity(PixelList,10);
    end
       
    %for small segments nan values might appear. set the value of nan entries to the mean of all non nan entries.

    nanInd = find(isnan(weights));
    [x,y] = ind2sub(size(weights),nanInd);
    for i = 1:length(x)
	weights(x(i),y(i)) = mean(weights(~isnan(weights(:,y(i))),y(i)));
    end

    weightNames = { 'shape: logObjSize', 'shape: pcFracs(1)', 'shape: pcFracs(2)', ...
        'shape: pcFracs(3)', 'shape: pcFracs(3)/pcFracs(2)', 'shape: pcFracs(3)/pcFracs(1)', ...
        'shape: pcFracs(2)/pcFracs(1)', 'shape: entropy' };
end

