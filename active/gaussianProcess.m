%% initialize covariance kernel
% squared exponential
kSigma = 1;
kLength =  1;
kFunc = @(x1,x2) kSigma^2*exp((x1-x2).^2/(-2*kLength^2));
% gaussian noise
nSigma = 0.1;
nFunc = @(x1,x2) nSigma^2*(x1==x2);
% sum together
k = @(x1,x2) kFunc(x1,x2)+nFunc(x1,x2);

%% look at gaussian process prior (?)
samplePoints=-1:0.01:1;
mean = zeros(1,length(samplePoints));
variance = k(samplePoints,samplePoints);
% plot using anonymus function
plot_variance = @(x,lower,upper,color) fill([x,x(end:-1:1)],[upper,fliplr(lower)],color,'EdgeColor',color);
plot_variance(samplePoints,mean-2*variance,mean+2*variance,[0.8 0.8 0.8]);
hold on;
plot(samplePoints,mean,'k','LineWidth',2);

%% sample from gaussian process
K_xx = k(repmat(samplePoints',1,length(samplePoints)),repmat(samplePoints,length(samplePoints),1));

[V,D]=eig(K_xx);
A=real(V*(D.^(1/2)));

for i=1:10
    gp_sample(:,i) = A * randn(length(samplePoints),1);
end

plot(samplePoints,real(gp_sample))

%% test gaussian process regression
X = [-1 -0.5 0 0.5 1];
Y = [2 1 0 1 2];

K_xx = k(repmat(X',1,length(X)),repmat(X,length(X),1));
K_ss = k(repmat(samplePoints',1,length(samplePoints)),repmat(samplePoints,length(samplePoints),1));
K_sx = k(repmat(samplePoints',1,length(X)),repmat(X,length(samplePoints),1));

mean = (K_sx/K_xx)*Y';
variance = 2*sqrt(diag(K_ss-K_sx/K_xx*K_sx'));

plot_variance(samplePoints,(mean-variance)',(mean+variance)',[0.8 0.8 0.8]);
hold on;
plot(samplePoints,mean,'k-','LineWidth',2);
plot(X,Y,'r.','MarkerSize',15);

evidence = exp((Y/K_xx*Y'+log(det(K_xx))+length(Y)*log(2*pi))/-2);
title (['evidence: ' num2str(evidence)]);

legend('confidence bounds','mean','data points','location','SouthEast');

%% sample from the gaussian process posterior
K_ss = k(repmat(samplePoints',1,length(samplePoints)),repmat(samplePoints,length(samplePoints),1));
K_sx = k(repmat(samplePoints',1,length(X)),repmat(X,length(samplePoints),1));

[V,D]=eig(K_ss-K_sx/K_xx*K_sx');
A=real(V*(D.^(1/2)));

for i=1:10
    gp_sample(:,i) = A * randn(length(A),1)+K_sx/K_xx*Y';
end
hold on;
plot(samplePoints,real(gp_sample));
plot(X,Y,'r.','MarkerSize',20);