function enforceMyelinSegments(p, newPrefix, bbox)
    % enforceMyelinSegments(p, box)
    %   The current membrane detection / segmentation pipeline
    %   struggles with myelin sheaths and produces 
    %   neuron-myelin mergers. This function uses the results of
    %   myelin detection and patches up the previously generated
    %   classification result to reduce the rate of mergers.
    %
    % p
    %   Parameter structure from pipeline repo for dataset
    %
    % bbox
    %   Bounding Box where the fix is to be applied
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    % Adpoted as first child by:
    %   Manuel
 
    % make sure that the bounding box
    % is aligned with the KNOSSOS cube hierarchy
    assert(Util.checkBoundingBox(bbox));
    
    % Load myelin detection
    myelin = loadSegDataGlobal(p.myelin, bbox);
    
    if isfield(p.myelin, 'segId')
        % NOTE(amotta): Manuel used to generate segmentations with the
        % heuristically detected blood vessels, nuclei, and myelin.
        myelin = myelin == p.myelin.segId;
    end
    
    myelin = logical(myelin);
    border = bwdist(~myelin, 'euclidean');
    border = (border > 0) & (border < 2);

    % Load classification data
    class = loadClassData(p.class, bbox);
    
    % Build up the wall
    borderHeight = -1.7;
    class(border) = borderHeight;
    
    % write result
    saveParam = Util.modifyStruct(p.class, 'prefix', newPrefix);
    saveOffset = reshape(bbox(:, 1), 1, []);
    saveClassData(saveParam, saveOffset, class);
end

