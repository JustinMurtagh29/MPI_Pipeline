function segmentation = watershedSeg_v2_retina( aff,cell)

tRange = cell{1};
vRange = cell{2};

N_LEVELS = 20000;
N_LEVELS_BOUNDARY = 1;

aff{1} = binAff(aff{1}, N_LEVELS);
aff{2} = binAff(aff{2}, N_LEVELS);
aff{3} = binAff(aff{3}, N_LEVELS);


segmentation = cell([length(tRange) length(vRange)]);

affAll = (aff{1} + aff{2} + aff{3}) / 3;

for t=1:length(tRange)
    display(num2str(t));

    for v=1:length(vRange)
    	marker = affAll < tRange(t);
        marker = bwareaopen(marker, vRange(v), 26);
        %segmentation{t,v} = watershed_3D(aff{1}, aff{2}, aff{3}, fgm1);
        
        marker = uint16(labelmatrix(bwconncomp(marker)));
        segmentation{t,v}=watershed_threeTimes3D(aff{1}, aff{2}, aff{3}, marker, N_LEVELS, N_LEVELS - N_LEVELS_BOUNDARY);
    end
end
clear affAll aff

end

