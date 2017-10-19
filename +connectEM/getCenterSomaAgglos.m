function somas = getCenterSomaAgglos

load('/gaba/u/mberning/results/pipeline/20170217_ROI/aggloState/somas_with_merged_somas.mat') 
somaCoordinates =   bsxfun(@times,[1392        2835        2779;
        3638        6240        1491;
        3044        1628         571;
        1573        2228        2355;
        2310        4859        1835;
        2098        6662         731;
        2278        2507         779;
        1237        6054        2563;
        2181        3208        2371;
        3797        5026        3035;
        3420        5468        2187;
        3363        6404        2867;
        3583        7443        2651;
        2162        1948        2659;
        1830        4736        2371;
        2727        7603         875;
        2801        5465        1595;
        2221        5974        2371;
        3064        2345        2571;
        3135        1758        1987;
        2559        4401        2859;
        4283        1483        1659;
        2384        1769        1723;
        3288        4610        2619;
        4422        3248        2075;
        1190        6023        1811;
        3966        3692        1331],[11.2400   11.2400   28]);
    
somaCOMs = bsxfun(@times,cell2mat(cellfun(@mean,somas(:,1),'uni',0)),[11.2400   11.2400   28]);
centerSomas = zeros(size(somaCoordinates,1),1);
for s = 1:size(somaCoordinates,1)
    [~,centerSomas(s)] = min(pdist2(somaCoordinates(s,:),somaCOMs));
end
assert(numel(unique(centerSomas))==size(somaCoordinates,1))
somas = somas(centerSomas,3);