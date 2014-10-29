w.Sigma = 1*eye(size(train.X,1));
[p, latMean, latVar, nlml] = binaryGPclassifier(w, train.X, train.y, train.X);

idx = train.y ~= 0;
a = p(idx);
b = p(~idx);

colors = [0 1 0; 1 0 0];
binSize = .01;
upperLimit = 1;
x = binSize/2:binSize:upperLimit-binSize/2;
c = hist(a, x)/numel(a);
d = hist(b, x)/numel(b);

figure;
h = bar(x, [c' d']);
for i=1:2
    set(h(i),'facecolor',colors(i,:),'edgecolor','none');
end

xlabel('probability of connectedness predicted by classifier');
legend({'connected' 'unconnected'});

