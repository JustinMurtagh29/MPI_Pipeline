function [estimator, log_like] = gaussianProbit( w, sigmaSquared, X, Y)

c_initial = (X'*X)\(X'*Y)*(rand-0.5)*10;
[estimator,log_like] = fminunc(@(c)ML_PROBIT(c,Y,X),c_initial);

function log_like = ML_PROBIT(c,Y,X)
q = 2*Y-1;
probit_F = normcdf(q.*(X*c));
log_like = sum(log(probit_F));
log_like = -log_like;
end

end
