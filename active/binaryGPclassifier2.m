function [out1, out2, out3, out4] = binaryGPclassifier2( prior, X, y, x )
% Ensure y = +/-1 if 0/1 is passed
y(y==0) = -1;
% Calculate kernel (everything still in data space so far)
K_XX = X'*prior.Sigma*X;
[a, sW, L, alml] = feval('laplaceApprox', K_XX, 'logistic', y);
if nargin==3
    % Mode-finding
    out1 = K_XX*a;
    out2 = alml;
    out3 = [];
    out4 = [];
else
    % Prediction
    nStar = size(x,1);
    nBatch = 10000;
    nBatchStart = 0;
    out1 = zeros(nStar,1);
    out2 = zeros(nStar,1);
    out3 = zeros(nStar,1);
    out4 = alml;
    while nBatchStart < nStar
        id = (nBatchStart+1):min(nBatchStart+nBatch,nStar);
        x_star = x(id,:)';
        K_xx = x_star'*prior.Sigma*x_star;
        K_Xx = X'*prior.Sigma*x_star;
        mu = K_Xx'*a;
        v  = L'\(repmat(sW,1,length(id)).*K_Xx);
        s2 = diag(K_xx) - sum(v.*v,1)';
        p  = feval('erfint', ones(size(mu)), mu, s2); % ones to find probabilities for class 1
        out1(id) = p;
        out2(id) = mu;
        out3(id) = s2;
        nBatchStart = id(end);
    end
end
end

function [a, sW, L, alml] = laplaceApprox(K, lik, y)
tol = 1e-9;
n = size(y,1);
f = zeros(n,1);
a = f;
% Calculate indicators for convergence criteria before loop (since it seems no until
% function is present in MATLAB)
[l,dl,d2l] = feval(lik,y,f); % Likelihood function & derivatives (logistic), see below
W = -d2l; % Hessian
psiNew = l;
psiOld = -Inf;
while psiNew - psiOld > tol
    psiOld = psiNew;
    aOld = a;
    sW = sqrt(W);
    L = chol(eye(n)+sW*sW'.*K); % L'*L=eye(n)+sW*K*sW
    b = W.*f+dl;
    a = b - sW.*(L\(L'\(sW.*(K*b))));
    f = K*a;
    [l,dl,d2l] = feval(lik,y,f);
    W = -d2l;
    psiNew = -a'*f/2 + l; % new likelihood
    % Suggestion in Rasmussen: Check whether objective decreases
    i = 0;
    while i < 10 && psiNew < psiOld
        a = (aOld+a)/2; % reduce step size by half otherwise
        f = K*a;
        [l,dl,d2l] = feval(lik,y,f);
        W = -d2l;
        psiNew = -a'*f/2 + l;
        i = i+1;
    end
    if psiNew < psiOld
        % If after 10 iterations still divergence, stop!
        error('Problem in newton iteration!');
    end
end
% Recalculate L, sW and alml to have correct output
sW = sqrt(W);
L  = chol(eye(n)+sW*sW'.*K);
alml = a'*f/2 - l + sum(log(diag(L))); % approximate log marginal likelihood
end

function [out1, out2, out3] = logistic(y, f)
% Function to calculate likelihood, log likelihood and derivatives for
% logistic (3 output arguments calculates derivates, otherwise
% likelihood & log likelihood)
	yf = y.*f;
    if nargout == 2
        lik = 1./(1+exp(-yf));
        logLik = log(lik);
        out1 = lik;
        out2 = logLik;
    else
        s    = -yf;
        ps   = max(0,s);
        out1 = -sum(ps+log(exp(-ps)+exp(s-ps)));
        nf    = min(0,f); 
        out2 = (y+1)/2-exp(nf)./(exp(nf)+exp(nf-f));
        out3 = -exp(2*nf-f)./(exp(nf)+exp(nf-f)).^2;
    end
end

function [m0,m1,m2] = erfint(y, mu, s2)
% Function to approximate Moments using 5 cumulative Gaussian distributions
% Taken from gpml toolbox by Rasmussen
l = [0.44 0.41 0.40 0.39 0.36]; % approximation coefficients lambda_i
c = [1.146480988574439e+02; -1.508871030070582e+03; 2.676085036831241e+03;  
    -1.356294962039222e+03;  7.543285642111850e+01                        ];
% zeroth moment    
S2 = 2*s2.*(y.^2)*(l.^2) + 1;                                    
S  = real(sqrt( S2 ));
Z  = mu.*y*l./S;
M0 = erf(Z);
m0 = ( 1 + M0*c )/2;   
% first moment
NormZ = exp(-Z.^2)/sqrt(2*pi);
M0mu = M0.*repmat(mu,[1,5]);
M1 = (2*sqrt(2)*y.*s2)*l.*NormZ./S + M0mu;
m1 = ( mu + M1*c )/2;
% second moment
M2 =   repmat(2*mu,[1,5]).*(1+s2.*y.^2*(l.^2)).*(M1-M0mu)./S2 ...
     + repmat(s2+mu.^2,[1,5]).*M0;
m2 = ( mu.^2 + s2 + M2*c )/2;
end
