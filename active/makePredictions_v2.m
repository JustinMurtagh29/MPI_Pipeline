function makePredictions_v2(w, X, Y, x_star, sigmaSquared, dim)
plot_variance = @(x,lower,upper,color) fill([x,x(end:-1:1)],[upper,fliplr(lower)],color,'EdgeColor',color);

x = zeros(size(w.Mu,1),size(x_star,2));
x(dim,:) = x_star;

K_XX = X'*w.Sigma*X + sigmaSquared * eye(size(X,2));
K_xx = x'*w.Sigma*x;
K_xX = x'*w.Sigma*X;

mean = (K_xX/K_XX)*Y';
variance = 2*sqrt(diag(K_xx-K_xX/K_XX*K_xX'));

plot_variance(x_star,(mean-variance)',(mean+variance)',[0.8 0.8 0.8]);
hold on;
plot(x_star,mean,'k-','LineWidth',2);
plot(X,Y,'r.','MarkerSize',15);

R = chol(K_XX); % R' * R = A; Y / (R'R) * Y = (Y/R) * (Y/R)'
YR = Y/R; halflogdetKXX = sum(log(diag(R)));
evidence = exp(-0.5*YR*YR'+2*halflogdetKXX+length(Y)*log(2*pi));
title (['evidence: ' num2str(evidence)]);

legend('confidence bounds','mean','data points','location','SouthEast');

end

