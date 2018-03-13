function [somas,somaCoords,KAMINsomaCoords] = getSomaAgglos(filename,type)
% returns the somata which overlap with the center, border or with all soma
% locations (which were extracted from the KAMIN list)
% INPUTS
% filename     string with path to one of the soma agglo mat files which
%              contains the Nx3 cell array "somas" (first column are the
%              agglo coordinates, third column are the agglo segIDs)
% type         string which can be center (default), border or all
%
% OUTPUTS
% somas        soma agglos in the old agglo format (cell array)
% somaCoords   all coordinates of each segId in each somaAgglo (scaled to nm)
% KAMINsomaCoords   the center coordinate from the KAMIN list (in pixels)
%
% by marcel.beining@brain.mpg.de

if ~exist('type','var') || isempty(type)
    type = 'center';
end
load(filename) 
centerSomaCoordinates =   [1392        2835        2779;
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
        3966        3692        1331];
    
    
allSomaCoordinates =   [1392	2835	2779;
128	1659	379;
3869	811	1171;
189	7128	2691;
5637	3960	1291;
3638	6240	1491;
5473	6700	51;
503	2404	779;
2258	3759	3235;
3044	1628	571;
1573	2228	2355;
2556	1933	3;
2321	8567	3051;
3673	3521	331;
164	6132	1811;
5162	7899	2811;
5294	6885	2075;
447	4609	2507;
3295	609	3411;
38	3888	451;
2440	2693	3155;
2310	4859	1835;
66	8146	1195;
63	5241	667;
3872	5306	395;
3980	7974	1075;
1734	527	1251;
2098	6662	731;
4502	5725	3411;
2278	2507	779;
4790	271	91;
1237	6054	2563;
2181	3208	2371;
2491	7932	1891;
510	2130	2371;
4998	5929	195;
5543	1548	1827;
4199	6135	363;
4923	7426	1011;
5176	3958	1867;
3797	5026	3035;
1389	1695	435;
3420	5468	2187;
1749	6991	91;
242	1157	2763;
1675	6899	3259;
1533	7733	2507;
3363	6404	2867;
4499	988	3091;
481	90	1859;
3525	8527	3411;
5545	6966	3307;
721	4059	1859;
5218	6918	2803;
416	7225	1803;
3583	7443	2651;
2162	1948	2659;
1830	4736	2371;
1457	7908	1523;
3591	484	3059;
5119	1971	2691;
5649	4365	1003;
2727	7603	875;
2801	5465	1595;
394	3508	955;
2539	4757	3379;
2221	5974	2371;
3064	2345	2571;
5024	2662	267;
3135	1758	1987;
5183	5679	1051;
153	3190	2291;
3410	69	2603;
636	370	3003;
2559	4401	2859;
1544	2477	3;
5507	1611	3171;
1412	8355	2675;
3412	6816	315;
4608	6937	3347;
157	2403	1251;
4283	1483	1659;
2384	1769	1723;
3832	2659	3155;
2460	2047	3387;
489	4956	1595;
3288	4610	2619;
2237	977	2371;
4422	3248	2075;
1190	6023	1811;
3966	3692	1331;
1778	565	739;
707	4030	499;
5199	4692	3339;
4678	4815	3411;
2761	7513	123;
3242	7673	243;
732	8436	1103;
5454	238	2519;
2099	7999	3397;
5497	3630	3332;
5482	450	3140;
239	7189	3366];

switch type
    case 'center'
        KAMINsomaCoords = centerSomaCoordinates;
    case 'all'
        KAMINsomaCoords = allSomaCoordinates;
    case 'border'
        KAMINsomaCoords = setdiff(allSomaCoordinates,centerSomaCoordinates,'rows');
    otherwise
        error('Type not found')
end
theseSomaCoordinatesScaled = bsxfun(@times, KAMINsomaCoords,[11.2400   11.2400   28]);

somaCOMs = bsxfun(@times,cell2mat(cellfun(@mean,somas(:,1),'uni',0)),[11.2400   11.2400   28]);
idxSomas = zeros(size(theseSomaCoordinatesScaled,1),1);
for s = 1:size(theseSomaCoordinatesScaled,1)
    [dst(s),idxSomas(s)] = min(pdist2(theseSomaCoordinatesScaled(s,:),somaCOMs));
end
if (numel(unique(idxSomas))~=size(theseSomaCoordinatesScaled,1))
    [~,ind,ind2] = unique(idxSomas);
    % get the KAMIN center coordinate which had the lowest distance to the
    % soma Agglo
    for n = 1:numel(ind)
        ind3 = find(n==ind2);
        [~,ind4] = min(dst(ind3));
        ind(n) = ind3(ind4);
    end
    warning('Caution! %d soma coordinates from the KAMIN list were part from the same soma agglo (probably because a somaAgglo does not exist for a KAMIN Id)!These are removed!',numel(idxSomas)-numel(ind))
    idxSomas = idxSomas(sort(ind));
    KAMINsomaCoords = KAMINsomaCoords(sort(ind),:);
end
somaCoords = cellfun(@(x) bsxfun(@times,x,[11.2400   11.2400   28]),somas(idxSomas,1),'uni',0);
somas = somas(idxSomas,3);
