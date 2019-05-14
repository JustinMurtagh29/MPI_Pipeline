agglo = load('/home/yy/segmentAggloPredictions-20170607.mat');

nonastro = agglo.segId(agglo.probs(:,1) <= 0.5);
astro = agglo.segId(agglo.probs(:,1) > 0.5);
components = { [0; nonastro], astro };

name = 'astroagglos-20170607';

folder = '/home/yy/';


WK.makeWKMapping( components, name , folder );