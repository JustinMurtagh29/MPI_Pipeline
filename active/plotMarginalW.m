function plotMarginalW( gaussian, dim, label )
% marginalize
mu = gaussian.Mu(dim);
sigma = gaussian.Sigma(dim,dim);

samplePoints = -3:0.01:3;
Y = normpdf(samplePoints,mu,sigma);

plot(samplePoints,Y,'k','LineWidth',2);
axis equal;
xlabel(['weight distribution for ' label{dim}]);

end

