function segIds = determineSegIdsFromMag4(p, mag4)
    % Extract segIds from globalized segmentation hit by label in KNOSSOS hierachy mag4

    % Load labels from Knossos hierachy and initalize
    labels = loadSegDataGlobal(mag4, mag4.bbox);
    segIds = cell(max(labels(:)), 1);
    cubesProcessed = 0;
    % Indices used for upsampling of labels to mag1 (start of if)
    x = ceil(0.25:0.25:(p.tileSize(1)/4));
    y = ceil(0.25:0.25:(p.tileSize(2)/4));
    z = ceil(0.25:0.25:(p.tileSize(3)/4));
    % Iterate over local segmentation cubes
    tic;
    for i=1:size(p.local,1)
        for j=1:size(p.local,2)
            for k=1:size(p.local,3)
                % Extrat label in local bounding box
                thisBBox = p.local(i,j,k).bboxSmall;
                thisBBox = thisBBox - repmat(p.bbox(:,1), 1,2) + ones(3,2);
                thisBBoxMag4 = ceil(thisBBox ./ 4);
                theseLabels = labels(thisBBoxMag4(1,1):thisBBoxMag4(1,2), ...
                    thisBBoxMag4(2,1):thisBBoxMag4(2,2), ...
                    thisBBoxMag4(3,1):thisBBoxMag4(3,2));
                if any(theseLabels(:))
                    % Upsample labels to mag1
                    theseLabels = theseLabels(x,y,z);
                    % Load bboxSmall part of segmentation
                    load([p.local(i,j,k).saveFolder 'segGlobal.mat']);
                    thisSeg = seg;
                    % Get linear indices of all labels
                    rp = regionprops(theseLabels(:), 'PixelIdxList');
                    % Find all segments for each label
                    for l=1:length(rp)
                        segIds{l} = unique([segIds{l}; thisSeg(rp(l).PixelIdxList)]);
                    end
                end
                cubesProcessed = cubesProcessed + 1;
                Util.progressBar(cubesProcessed, prod(size(p.local)));
            end
        end
    end

    % Remove zeroes from agglomeration
    segIds = cellfun(@(x)x(x~=0), segIds, 'UniformOutput', false);

end

