function enforceMyelinMasks(p, newClassRoot, bbox)
    % enforceMyelinSegments(p, box)
    %   The current membrane detection / segmentation pipeline
    %   struggles with myelin sheaths and produces 
    %   neuron-myelin mergers. This function uses the results of
    %   myelinCNN detection and patches up the previously generated
    %   classification result to reduce the rate of mergers.
    %
    % p
    %   Parameter structure from pipeline repo for dataset
    %
    % bbox
    %   Bounding Box where the fix is to be applied
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    %   Sahil Loomba <sahil.loomba@brain.mpg.de>
    % Adpoted as first child by:
    %   Manuel
 
    % make sure that the bounding box
    % is aligned with the KNOSSOS cube hierarchy
    assert(Util.checkBoundingBox(bbox));
    
    % Load myelin detection
    myelin = loadClassData(p.classMyelin, bbox);
   
    % Threshold myelin voxels at thr 
    myelin = myelin > p.myelin.thr; 
    border = bwdist(~myelin, 'euclidean');
    border = (border > 0) & (border < 2);

    % Load classification data
    class = loadClassData(p.class, bbox);
    
    % Build up the wall
    borderHeight = -1.7;
    class(border) = borderHeight;
    
    % write result
    saveParam = Util.modifyStruct(p.class, 'root', newClassRoot, 'backend', 'wkwrap');
    saveOffset = reshape(bbox(:, 1), 1, []);
    saveClassData(saveParam, saveOffset, class);
end

