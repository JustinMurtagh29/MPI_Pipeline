function segmentation = watershedSeg_v3_cortex( aff, cell)
%v3 watershedSeg with AnisoDiff 
%  kRange=[], niter=[0:100], option=[1,2,3]                                                                                                
diff = single(aff);

hRange = cell{1};
vRange = cell{2};
kRange = cell{3};
niter = cell {4};
option = cell {5};

segmentation = cell([length(hRange) length(vRange) length(kRange) length(niter) length(option)]);

for h=1:length(hRange)
    for v=1:length(vRange)
        for k=1:length(kRange)
            for i=1:length(niter)
                for o=1:length(option)
                    
                diff = AnisoDiff( diff, kRange, niter, option);        
                affHmin= imhmin(diff, hRange(h), 26);
                bw1 = imregionalmin(affHmin, 26); 
                
                end
            end
        end
    end 
    bw1 = bwareaopen(bw1, vRange(v), 26);
    affImposed = imimposemin(aff, bw1);
    segmentation{h,v} = watershed(affImposed, 26);      
end

end