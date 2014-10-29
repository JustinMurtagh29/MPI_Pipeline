function plotMarginal( gaussian, x, xY, Y, label )
% pass gaussian struct and two dimensions
plot_variance = @(x,lower,upper,color) fill([x,x(end:-1:1)],[upper,fliplr(lower)],color,'EdgeColor',color);

mean = gaussian.Mu';
variance = diag(gaussian.Sigma)';

plot_variance(x,mean-variance,mean+variance,[0.8 0.8 0.8 ]);
hold on;
plot(x,gaussian.Mu,'k','LineWidth',2);
plot(xY,Y,'r.','MarkerSize',15);
axis equal;
xlabel(['weight distribution for ' label]);

end