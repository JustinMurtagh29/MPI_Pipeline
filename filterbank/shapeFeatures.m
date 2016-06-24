function [weights, weightNames] = shapeFeatures(segments,siz)

    for i=1:length(segments)
        % switch from ind to sub
        [x,y,z] = ind2sub(siz,double(segments{i}));
        PixelList = [x y z];             
        objSize = size(PixelList,1);
        % added MB, clean up, getting memory errors on fermat
        clear x y z;

        % Hmm not sure why Benjamin chose PC (and only latent output), extend?
        [COEFF,SCORE,latent] = princomp(PixelList, 'econ');
        pcFracs = latent./sum(latent);
        
        if objSize > 3
            weights(i,1) = log(objSize);                     
            weights(i,2) = pcFracs(1);
            weights(i,3) = pcFracs(2);
            weights(i,4) = pcFracs(3);
            weights(i,5) = pcFracs(3)/pcFracs(2);
            weights(i,6) = pcFracs(3)/pcFracs(1);
            weights(i,7) = pcFracs(2)/pcFracs(1);
            weights(i,8) = -pcFracs(1)*log(pcFracs(1))-pcFracs(2)*log(pcFracs(2))-pcFracs(3)*log(pcFracs(3)); % entropy
        else
            weights(i,:) = 0;
        end
       
    %for small segments nan values might appear. set the value of nan entries to zero 
    nanInd = find(isnan(weights));
    [x,y] = ind2sub(size(weights),nanInd);
    for i = 1:length(x)
	    weights(x(i),y(i)) = 0;
    end

    weightNames = { 'shape: logObjSize', 'shape: pcFracs(1)', 'shape: pcFracs(2)', ...
        'shape: pcFracs(3)', 'shape: pcFracs(3)/pcFracs(2)', 'shape: pcFracs(3)/pcFracs(1)', ...
        'shape: pcFracs(2)/pcFracs(1)', 'shape: entropy' };
end

