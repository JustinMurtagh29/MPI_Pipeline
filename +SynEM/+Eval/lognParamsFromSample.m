function [ mu, sig ] = lognParamsFromSample( m, v )
%lognParamsFromSample Convert mean and variance of a non-logarithmized
% sample to the lognormal distribution parameters mu and sig.
% INPUT m: double
%           Mean of the sample.
%       v: double
%           Variance of the sample.
% OUTPUT mu: double
%           Mu parameter of the lognormal distribution.
%        sig: double
%           Sigma parameter of the lognormal distribution.
%
% EXAMPLE
% mu = 3.2; sig = 1.6;
% a = lognrnd(mu, sig, 1e4, 1);
% [mu2, sig2] = lognParamsFromSample(mean(a), var(a));
% %mu2 and sig2 are the estimated for mu and sig
%
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

mu = log(m^2/sqrt(v + m^2));
sig = sqrt(log(v/m^2 + 1));


end

