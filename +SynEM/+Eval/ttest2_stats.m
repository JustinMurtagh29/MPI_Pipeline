function [ h, p ] = ttest2_stats( x, mu_y, s_y, m )
%TTEST2_stats Two sample two-sided t-test where the second distribution is
% given by its summary statistics.
% INPUT x: [Nx1] float
%           Samples of the first distribution.
%       mu_y: float
%           Mean of the second distribution.
%       s_y: float
%           Standard deviation of the second distribution.
%       m: int
%           Sample size of the second distibution.
% OUTPUT h: logical
%           Test decision for the rejection of the null hypothesis that x
%           and y come from independent random sample from a normal
%           distribution with equal mean and equal but unknown variance.
%           1 means that the null hypothesis is rejected.
%        p: float
%           p value of the test.
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

alpha = 0.05;
n = length(x);
mu_x = mean(x);
s_x = std(x);
df = n + m - 2;
s = sqrt(((n-1)*s_x^2 + (m-1)*s_y^2) / df);
t = sqrt(n*m/(n+m))*(mu_x - mu_y)/s;
p = 2 * tcdf(-abs(t), df);
h = p < alpha;

end

