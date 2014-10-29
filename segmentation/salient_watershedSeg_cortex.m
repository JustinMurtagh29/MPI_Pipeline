function segmentation = salient_watershedSeg_cortex(aff, cell)
%salient watershed with canny/Pb Filter

% aff = load('P:\20130603cube\aff\20130603classifierNew.mat');
% diff = single(imcomplement(aff.class(1:640,1:768,1:2)));

hRange = cell{1};
vRange = cell{2};
sigma = cell {3};
th_low = cell {4};
th_up = cell {5};

segmentation = cell([length(hRange) length(vRange) length(sigma) length(th_low) length(th_up)]);

for h=1:length(hRange)
    for v=1:length(vRange)
        for s=1:length(sigma)
            for t_l = 1:length(th_low)
                for t_u = 1:length(th_up)
                    
%                     mean = ordfilt3D(aff, 14, 'both');
                    canny = imcomplement(canny3D(aff, 3, sigma(s), th_up(t_u),th_low(t_l)));     
                    Pb = Pbfilter(aff); clear mean;

                    gradSal = canny.*Pb; clear canny Pb;
                %     [sal, IDX] = bwdist(gradSal,'euclidean');
                    sal = exp(-2.*gradSal); clear gradSal

                    bw1 = bwareaopen(imregionalmax(imhmin(sal, hRange(h), 26), 26), vRange(v), 26);
                    affImposed = imimposemin(aff, bw1);

                    segmentation{h,v,s,t_l,t_u} = watershed(affImposed,26);
    
                end
            end
        end
    end
end

