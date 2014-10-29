function plotMarginal2D( gaussian, dims, label )
% marginalize
mu = gaussian.Mu(dims);
sigma = gaussian.Sigma(dims,dims);
% eigendecomposition
[V,D]=eig(sigma);
A=V*(D.^(1/2)); A = chol(sigma + 1.0e-9 * eye(size(sigma)))';
% sample
samples = zeros(2,1000);
for i=1:1000
    samples(:,i) = A * randn(2,1) + mu;
end
% plot
plot(samples(1,:),samples(2,:),'.');
hold on;
plot(mu(1), mu(2),'xr','MarkerSize',10,'LineWidth',2);
axis equal;
xlabel(label(dims(1)));
ylabel(label(dims(2)));
end

